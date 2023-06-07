# python dependencies
from pathlib import Path
import gzip
from mimetypes import guess_type
from functools import partial
from uuid import uuid4
import zlib

# external dependencies
from Bio import SeqIO, Seq
from io import StringIO

# internal dependencies
from socialgene.utils.logging import log
import socialgene.utils.file_handling as fh
import re
from collections import Counter

# see SequenceParser at bottom for main class


def get_seqio_start(seq_feature):
    return seq_feature.location.start.real + 1


def get_seqio_end(seq_feature):
    return seq_feature.location.end.real


class GenbankParser:
    def __init__(self):
        self._self._count_loci_in_file
        pass

    def _add_assembly(self, input_path, seq_record):
        if hasattr(seq_record, "dbxrefs"):
            # try to retrieve assembly id from the seqrecord dbxref
            assembly_id = self._extract_from_dbxrefs(
                input_list=getattr(seq_record, "dbxrefs", None), prefix="Assembly:"
            )
        if not assembly_id:
            # use filename as assembly id
            assembly_id = str(Path(str(input_path).removesuffix(".gz")).stem)
        if not assembly_id:
            # use random unique id as assembly id
            assembly_id = str(uuid4())
        self.add_assembly(id=assembly_id)
        return assembly_id

    def _add_locus(self, seq_record, assembly_id):
        # add a locus object to an assembly object
        self.assemblies[assembly_id].add_locus(id=seq_record.id)
        return seq_record.id

    def _add_taxon(self, seq_feature, assembly_id):
        try:
            self.assemblies[assembly_id].taxid = self._extract_from_dbxrefs(
                input_list=seq_feature.qualifiers.get("db_xref"), prefix="taxon:"
            )
            log.debug(f"{assembly_id} is taxon: {self.assemblies[assembly_id].taxid}")
        except Exception as e:
            log.info(e)

    def _process_feature_note(self, seq_feature):
        note = None
        if "note" in seq_feature.qualifiers:
            note = seq_feature.qualifiers["note"]
            if isinstance(note, list):
                note = note[0]
            # Remove long note about gene prediction
            to_replace = "Derived by automated computational analysis using gene prediction method.*"
            note = re.sub(to_replace, "", note)
            if note == "":
                note = None
        return note

    def _process_go(self, seq_feature):
        gos = {"GO_component": None, "GO_process": None, "GO_function": None}
        for k in gos.keys():
            if k in seq_feature.qualifiers:
                # e.g. -> seq_feature.qualifiers[k]=['GO:0006260 - DNA replication [Evidence IEA]']
                if k in seq_feature.qualifiers:
                    gos[k] = ";".join(
                        [i.split(" ", 1)[0] for i in seq_feature.qualifiers[k]]
                    )
        return {k.lower(): v for k, v in gos.items()}

    def _add_feature(
        self,
        seq_record,
        seq_feature,
        assembly_id,
        locus_id,
        keep_sequence,
    ):
        # Grab the locus/protein description if it exists
        if "product" in seq_feature.qualifiers:
            product = seq_feature.qualifiers["product"][0]
        else:
            product = None
        if seq_feature.type == "source":
            # TODO: should this overwrite (current) or compare every iteration per assembly?
            for k, v in seq_feature.qualifiers.items():
                if k in self.assemblies[assembly_id].info:
                    self.assemblies[assembly_id].info[k] = v
                    self.assemblies[assembly_id].loci[locus_id].info[k] = v
        if seq_feature.type in ["protein", "CDS"]:
            try:
                # Assign a protein id if it exists, a locus tag if it doesn't (e.g. pseudogene)
                # if neither, make something up
                locus_tag = None
                if "protein_id" in seq_feature.qualifiers:
                    protein_id = seq_feature.qualifiers["protein_id"]
                elif "locus_tag" in seq_feature.qualifiers:
                    # use locus tag as protein instead
                    protein_id = seq_feature.qualifiers["locus_tag"]
                    locus_tag = seq_feature.qualifiers["locus_tag"][0]
                    if isinstance(protein_id, list):
                        locus_tag = locus_tag[0]
                else:
                    protein_id = uuid4()

                if isinstance(protein_id, list):
                    protein_id = protein_id[0]
                if "translation" in seq_feature.qualifiers:
                    translation = seq_feature.qualifiers["translation"][0]
                elif "pseudo" or "pseudogene" in seq_feature.qualifiers:
                    # removes biopython potential deprecation warning to prepend Ns if non-modulo 3 seqs,
                    # may be brittle if biopython changes
                    translation = str(
                        Seq(seq_feature.extract(seq_record).seq).translate()
                    )
                    product = f"pseudo_{product}"
                    protein_id = seq_feature.qualifiers["locus_tag"][0]
                else:
                    raise ValueError(
                        f"Panic!!! Not a protein or pseudo protein: {seq_feature.qualifiers}"
                    )
                description = product.strip()
                hash_id = self.add_protein(
                    description=description,
                    external_protein_id=protein_id.strip(),
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
                    locus_tag=locus_tag,
                    description=description,
                    note=self._process_feature_note(seq_feature),
                    **self._process_go(seq_feature),
                )
            except Exception as e:
                log.debug(e)

    def _parse_record(self, seq_record, assembly_id, keep_sequence):
        locus_id = self._add_locus(seq_record=seq_record, assembly_id=assembly_id)
        for seq_feature in seq_record.features:
            self._add_feature(
                seq_record=seq_record,
                seq_feature=seq_feature,
                assembly_id=assembly_id,
                locus_id=locus_id,
                keep_sequence=keep_sequence,
            )
        self._count_loci_in_file = self._count_loci_in_file + Counter(
            [i.type for i in seq_record.features]
        )

    def _parse_genbank(
        self,
        handle,
        input_path,
        assembly_id=None,
        keep_sequence=True,
    ):
        """Internal function to parse a read genbank file, per SEQIO record
        Args:
            handle (StringIO): genbank file handle
            input_path (str): Path to genbank file (filename used as assembly name, if missing)
            assembly_id (str, optional): Assign the assembly id (used for for neo4j)
            keep_sequence (bool, optional): Store the amino acid sequence of proteins?
        """
        # only set assembly info on reading first record
        stop_after_one = True
        for seq_record in SeqIO.parse(handle, "genbank"):
            # Add assembly
            if stop_after_one:
                if not assembly_id:
                    assembly_id = self._add_assembly(
                        input_path=input_path, seq_record=seq_record
                    )
                self._add_taxon(
                    seq_feature=seq_record.features[0], assembly_id=assembly_id
                )
                stop_after_one = False
            # Add locus and its features
            log.debug(f"Processing assembly: {assembly_id}")
            self._parse_record(
                seq_record=seq_record,
                assembly_id=assembly_id,
                keep_sequence=keep_sequence,
            )
        log.debug(f"Processing seq_record: {seq_record.id}")

    def _open_genbank(self, input_path, **kwargs):
        """Parse a genbank file
        Args:
            input_path (str): path to genbank file (should be able to handle .gz and tar archives)
        """
        # self._count_loci_in_file is mostly for sending a summary of parsed features to the logger
        self._count_loci_in_file = Counter()
        input_path = Path(input_path)
        if fh.check_if_tar(input_path):
            tar = fh.opener_result(input_path)
            log.info(f"tar archive has {len(tar.members)} files")
            log.info(tar)
            for member in tar:
                f = tar.extractfile(member)
                if "gbk" in member.name or "gbff" in member.name:
                    try:
                        log.info(f"Extracting and parsing: {member.name}")
                        input_ = zlib.decompress(f.read(), 16 + zlib.MAX_WBITS).decode(
                            "utf-8"
                        )
                        with StringIO(input_) as handle:
                            self._parse_genbank(
                                handle=handle,
                                input_path=input_path,
                            )
                    except Exception as e:
                        log.debug(e)
        else:
            with fh.open_file(input_path) as handle:
                self._parse_genbank(
                    handle=handle,
                    input_path=input_path,
                    **kwargs,
                )
        log.info(
            f"'{input_path}' features {dict(sorted(dict(self._count_loci_in_file).items(), key=lambda item: item[1], reverse=True))}"
        )

    @staticmethod
    def _extract_from_dbxrefs(input_list, prefix):
        """Helper function for genbank parser, extract a prefixed string from list

        Args:
            input_list (list): e.g. ['BioProject:PRJNA224116', 'BioSample:SAMN00006196', 'Assembly:GCF_000145235.1']
            prefix (str): e.g. "BioProject:"

        Returns:
            list: "PRJNA224116"
        """
        temp = [x.removeprefix(prefix) for x in input_list if x.startswith(prefix)]
        return ";".join(temp)


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
                    external_protein_id=seq_record.id,
                    sequence=str(seq_record.seq),
                )
                self.assemblies[assembly_id].loci[assembly_id].add_feature(
                    type="protein", id=hash_id, start=0, end=0, strand=0
                )
                record_counter += 1
                count_proteins_in_file += 1
        log.info(f"Read {dict(count_proteins_in_file)} proteins from {input}")

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
                    external_protein_id=seq_record.id,
                    sequence=str(seq_record.seq),
                )
                self.assemblies[assembly_id].loci[assembly_id].add_feature(
                    type="CDS",
                    id=prot_hash,
                )
                record_counter += 1
                count_proteins_in_file += 1
        log.info(f"Read {dict(count_proteins_in_file)} proteins from input string")


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
            self._open_genbank(filepath, **kwargs)
        elif filetype == "fasta":
            self.parse_fasta_file(filepath)
        else:
            raise NotImplementedError(
                "May not be implemented, or you need to use the genbank/fasta parser directly. (e.g. for tar archives)"
            )
