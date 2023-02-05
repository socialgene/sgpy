# python dependencies
from pathlib import Path
import gzip
from mimetypes import guess_type
from functools import partial
from uuid import uuid4
import zlib

# external dependencies
from Bio import SeqIO
from io import StringIO


# internal dependencies
from socialgene.utils.logging import log
import socialgene.utils.file_handling as fh

# see SequenceParser at bottom for main class


def get_seqio_start(seq_feature):
    return seq_feature.location.start.real + 1


def get_seqio_end(seq_feature):
    return seq_feature.location.end.real


class GenbankParser:
    # TODO: https://github.com/kblin/ncbi-acc-download/blob/master/ncbi_acc_download/validate.py
    def __init__(self):
        pass

    def _parse_genbank(
        self,
        input_path,
        seq_record,
        count_loci_in_file={},
        total_locus_counter=0,
        assembly_id=None,
        keep_sequence=True,
    ):
        """Internal function to parse a read genbank file, per SEQIO record

        Args:
            input_path (str): Path to genbank file
            seq_record (SeqRecord): Biopython SeqRecord
            count_loci_in_file (dict, optional): Defaults to {}.
            total_locus_counter (int, optional): Defaults to 0.
            assembly_id (str, optional): Assign the assembly id (used for for neo4j)
        """
        if not assembly_id:
            try:
                # try to retrieve assembly id from the seqrecord object
                assembly_id = [
                    i for i in seq_record.dbxrefs if i.startswith("Assembly:")
                ][0]
                assembly_id = assembly_id.replace("Assembly:", "", 1)
                self.add_assembly(id=assembly_id)
            except Exception:
                try:
                    assembly_id = str(Path(str(input_path).removesuffix(".gz")).stem)
                except Exception:
                    assembly_id = str(uuid4())
                self.add_assembly(id=assembly_id)
        # because locus ids aren't guaranteed to be unique, esp across loci /files
        # track_locus_ids = {}
        for seq_feature in seq_record.features:  # feature
            # Add locus
            locus_id = seq_record.id
            self.assemblies[assembly_id].add_locus(id=locus_id)
            # keep track of the number of different features
            if seq_feature.type in count_loci_in_file:
                count_loci_in_file[seq_feature.type] += 1
            else:
                count_loci_in_file[seq_feature.type] = 1
            # Assign a protein id if it exists, a locus tag if it doesn't (e.g. pseudogene)
            # if neither, make something up
            if "protein_id" in seq_feature.qualifiers:
                protein_id = seq_feature.qualifiers["protein_id"]
            elif "locus_tag" in seq_feature.qualifiers:
                protein_id = seq_feature.qualifiers["locus_tag"]
            else:
                protein_id = uuid4()
            if isinstance(protein_id, list):
                protein_id = protein_id[0]
            # Grab the locus/protein description if it exists
            if "product" in seq_feature.qualifiers:
                product = seq_feature.qualifiers["product"][0]
            else:
                product = None
            if seq_feature.type == "source":
                try:
                    taxon_list = [
                        i
                        for i in seq_feature.qualifiers.get("db_xref")
                        if i.startswith("taxon:")
                    ]
                    taxon_str = taxon_list[0]
                    taxon_str = taxon_str.replace("taxon:", "")
                    self.assemblies[assembly_id].taxid = int(taxon_str)
                except Exception:
                    pass
                # TODO: should this overwrite (current) or compare every iteration per assembly?
                for k in self.assemblies[assembly_id].info.keys():
                    if k in seq_feature.qualifiers:
                        self.assemblies[assembly_id].info[k] = seq_feature.qualifiers[k]
                        self.assemblies[assembly_id].loci[locus_id].info[
                            k
                        ] = seq_feature.qualifiers[k]
            if any([True for i in ["protein", "CDS"] if i == seq_feature.type]):
                try:
                    if "translation" in seq_feature.qualifiers:
                        translation = seq_feature.qualifiers["translation"][0]
                    elif "pseudo" or "pseudogene" in seq_feature.qualifiers:
                        translation = str(
                            seq_feature.extract(seq_record).seq.translate()
                        )
                        product = f"pseudo_{product}"
                        protein_id = seq_feature.qualifiers["locus_tag"][0]
                    else:
                        raise ValueError(
                            f"Panic!!! Not a protein or pseudo protein: {seq_feature.qualifiers}"
                        )
                    hash_id = self.add_protein(
                        description=product.strip(),
                        other_id=protein_id.strip(),
                        sequence=translation.strip(),
                    )
                    if not keep_sequence:
                        self.proteins[hash_id].sequence = None

                    self.assemblies[assembly_id].loci[locus_id].add_feature(
                        type=seq_feature.type,
                        id=hash_id,
                        start=get_seqio_start(seq_feature),
                        end=get_seqio_end(seq_feature),
                        strand=seq_feature.location.strand,
                    )
                except Exception as e:
                    _ = e
                    pass
            else:
                # use incremented counter to id non-protein loci
                total_locus_counter += 1

    def _read_genbank(self, input_path, **kwargs):
        """Parse a genbank file
        Args:
            input_path (str): path to genbank file (should be able to handle .gz and tar archives)
        """
        count_loci_in_file = {}
        total_locus_counter = 0
        input_path = Path(input_path)
        # note: seq_record is LOCUS in a genbank file
        if fh.check_if_tar(input_path):
            tar = fh.opener_result(input_path)
            log.info(f"tar archive has {len(tar.members)} files")
            log.info(tar)
            for member in tar:
                f = tar.extractfile(member)
                # TODO: fix reliance on extension
                if member.name.endswith(".gbk.gz") or member.name.endswith(".gbff.gz"):
                    log.info(f"Extracting and parsing: {member.name}")
                    input_ = zlib.decompress(f.read(), 16 + zlib.MAX_WBITS).decode(
                        "utf-8"
                    )
                    with StringIO(input_) as handle:
                        for seq_record in SeqIO.parse(handle, "genbank"):
                            self._parse_genbank(
                                input_path=input_path,
                                seq_record=seq_record,
                                count_loci_in_file=count_loci_in_file,
                                total_locus_counter=total_locus_counter,
                            )
                else:
                    try:
                        input_ = f.read().decode("utf-8")
                        log.info(f"Extracting: {member.name}")
                        with StringIO(input_) as handle:
                            for seq_record in SeqIO.parse(handle, "genbank"):
                                self._parse_genbank(
                                    input_path=input_path,
                                    seq_record=seq_record,
                                    count_loci_in_file=count_loci_in_file,
                                    total_locus_counter=total_locus_counter,
                                )
                    except Exception:
                        pass
        else:
            with fh.open_file(input_path) as handle:
                for seq_record in SeqIO.parse(handle, "genbank"):
                    self._parse_genbank(
                        input_path=input_path,
                        seq_record=seq_record,
                        count_loci_in_file=count_loci_in_file,
                        total_locus_counter=total_locus_counter,
                        **kwargs,
                    )
        log.info(f"Read {count_loci_in_file} from {input_path}")

    @staticmethod
    def extract_from_dbxrefs(input_list, prefix):
        """Helper function for genbank parser, extract a prefixed string from list

        Args:
            input_list (list): e.g. ['BioProject:PRJNA224116', 'BioSample:SAMN00006196', 'Assembly:GCF_000145235.1']
            prefix (str): e.g. "BioProject:"

        Returns:
            list: "PRJNA224116"
        """
        temp = [x.removeprefix(prefix) for x in input_list if x.startswith(prefix)]
        if len(temp) == 0:
            temp = None
        else:
            temp = temp[0]
        return temp


class FastaParser:
    def __init__(self):
        pass

    # havent tested/worked on this since lots of updates
    def parse_fasta_file(self, input=None):
        """Parse a protein fasta file

        Args:
            input (str): path to fasta file
        """
        # if not appending, reset self.records
        record_counter = 0
        count_proteins_in_file = 0
        input = Path(input)
        try:
            assembly_id = input.stem
        except Exception:
            assembly_id = str(uuid4())
        encoding = guess_type(input)[1]
        if encoding == "gzip":
            _open = partial(gzip.open, mode="rt")
        else:
            _open = open
        with _open(input) as handle:
            self.add_assembly(id=assembly_id)
            self.assemblies[assembly_id].add_locus(id=assembly_id)
            for seq_record in SeqIO.parse(handle, "fasta"):
                # probably could use some more defensive programming here
                hash_id = self.add_protein(
                    description=seq_record.description,
                    other_id=seq_record.id,
                    sequence=str(seq_record.seq),
                )
                self.assemblies[assembly_id].loci[assembly_id].add_feature(
                    type="protein", id=hash_id, start=0, end=0, strand=0
                )
                record_counter += 1
                count_proteins_in_file += 1
        log.info(f"Read {count_proteins_in_file} proteins from {input}")

    def parse_fasta_string(self, input):
        """Parse a protein fasta file from a string

        Args:
            input (str): path to fasta file
        """
        # if not appending, reset self.records
        record_counter = 0
        count_proteins_in_file = 0
        with StringIO(input) as handle:
            assembly_id = str(Path(input).stem)
            self.add_assembly(id=assembly_id)
            self.assemblies[assembly_id].add_locus(id=assembly_id)
            for seq_record in SeqIO.parse(handle, "fasta"):
                prot_hash = self.add_protein(
                    description=seq_record.description,
                    other_id=seq_record.id,
                    sequence=str(seq_record.seq),
                )
                self.assemblies[assembly_id].loci[assembly_id].add_feature(
                    type="CDS",
                    id=prot_hash,
                    start=0,
                    end=0,
                    strand=0,
                )
                record_counter += 1
                count_proteins_in_file += 1
        log.info(f"Read {count_proteins_in_file} proteins from input string")


class SequenceParser(GenbankParser, FastaParser):
    def __init__(self):
        super().__init__()

    def parse(self, filepath, **kwargs):
        """Parse sequence files (main function)

        Args:
            filepath (str): path to sequence file
        """
        # TODO: fh.check_if_tar(filepath=filepath)
        filetype = fh.guess_filetype(filepath)
        if filetype == "genbank":
            self._read_genbank(filepath, **kwargs)
        elif filetype == "fasta":
            self.parse_fasta_file(filepath)
        else:
            raise NotImplementedError(
                "May not be implemented, or you need to use the genbank/fasta parser directly. (e.g. for tar archives)"
            )
