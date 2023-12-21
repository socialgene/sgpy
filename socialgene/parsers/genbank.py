import re
import zlib
from collections import Counter
from io import StringIO
from pathlib import Path
from typing import Dict, List
from uuid import uuid4

from Bio import Seq, SeqIO

import socialgene.utils.file_handling as fh
from socialgene.utils.logging import log


def get_seqio_start(seq_feature):
    return seq_feature.location.start.real + 1


def get_seqio_end(seq_feature):
    return seq_feature.location.end.real


class GenbankParserMixin:
    def __init__(self):
        self._self._count_loci_in_file
        pass

    def _add_assembly(self, input_path, seq_record):
        """
        The function `_add_assembly` adds an assembly ID to a sequence record, using either the dbxrefs
        attribute, the input file name, or a random unique ID.

        Args:
          input_path: The input_path parameter is the path to the input file. It is used to determine
        the assembly ID if it is not provided in the seq_record object.
          seq_record: The `seq_record` parameter is an object that represents a sequence record. It
        likely contains information about a biological sequence, such as DNA or protein, including its
        sequence data and associated metadata.

        Returns:
          the assembly ID.
        """
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
        self.add_assembly(uid=assembly_id, parent=self)
        return assembly_id

    def _add_locus(self, seq_record, assembly_id):
        """
        The function adds a locus object to an assembly object and returns the ID of the added sequence
        record.

        Args:
          seq_record: The `seq_record` parameter is a sequence record object. It typically contains
        information about a biological sequence, such as its ID, sequence data, and annotations. In this
        context, it is being used to add a locus object to an assembly object.
          assembly_id: The `assembly_id` parameter is the identifier of the assembly object to which the
        locus object will be added.

        Returns:
          the ID of the sequence record that was added to the assembly object.
        """
        # add a locus object to an assembly object
        self.assemblies[assembly_id].add_locus(external_id=seq_record.id)
        return seq_record.id

    def _add_taxon(self, seq_feature, assembly_id):
        """
        The function `_add_taxon` extracts the taxon ID from a sequence feature's qualifiers and assigns
        it to the corresponding assembly ID in a dictionary.

        Args:
          seq_feature: The `seq_feature` parameter is a sequence feature object that represents a
        specific feature of a DNA or protein sequence, such as a gene or a coding region.
          assembly_id: The `assembly_id` parameter is an identifier for a specific assembly. It is used
        to access and modify the `taxid` attribute of the assembly object in the `self.assemblies`
        dictionary.
        """
        try:
            self.assemblies[assembly_id].taxid = self._extract_from_dbxrefs(
                input_list=seq_feature.qualifiers.get("db_xref"), prefix="taxon:"
            )
            log.debug(f"{assembly_id} is taxon: {self.assemblies[assembly_id].taxid}")
        except Exception as e:
            log.debug(e)

    def _process_feature_note_protein_exceptions(self, note: str) -> Dict[str, bool]:
        """
        Processes a given feature note and returns a dictionary indicating whether certain protein exceptions are present.

        Args:
          note (str): The `note` parameter is a string that contains information about a locus feature.

        Returns:
          a dictionary where the keys are the keys from the `bad_proteins` dictionary and the values are
        the string "true" for each key that matches a pattern in the `note` string.
        """
        bad_proteins = {
            "partial_on_complete_genome": "partial on complete genome",
            "missing_start": "missing start",
            "missing_stop": "missing stop",
            "internal_stop": "internal stop",
            "partial_in_the_middle_of_a_contig": "partial in the middle of a contig",
            "missing_N_terminus": "missing N-terminus",
            "missing_C_terminus": "missing C-terminus",
            "frameshifted": "frameshifted",
            "too_short_partial_abutting_assembly_gap": "too short partial abutting assembly gap",
            "incomplete": "incomplete",
        }
        # The exact string "true" is looked for by neo4j admin import, not True
        # https://neo4j.com/docs/operations-manual/current/tools/neo4j-admin/neo4j-admin-import/#import-tool-header-format-properties
        return {k: "true" for k, v in bad_proteins.items() if re.search(v, note)}

    def _process_go(self, note: str) -> Dict[str, List[str]]:
        """
        Processes a given feature note and returns a dictionary with a single
        key "goterms" whose value is a list of strings representing GO terms found in the input string.

        Args:
          note (str): The `note` parameter is a string that contains information about a locus feature.

        Returns:
          a dictionary with a single key "goterms" and a value that is a list of strings. The strings in
        the list are extracted from the input string "note" using regular expression pattern matching.
        The regular expression pattern "GO:[0-9]{7}" matches strings that start with "GO:" followed by
        exactly 7 digits.
        """
        return {"goterms": re.findall("GO:[0-9]{7}", note)}

    def _process_feature_note(self, seq_feature) -> Dict:
        """
        The function `_process_feature_note` processes the note qualifiers of a sequence feature and
        returns a dictionary containing information extracted from the note.

        Args:
          seq_feature: The `seq_feature` parameter is a SeqFeature object that represents a feature
        (e.g., gene, CDS, protein) in a sequence record. It contains information such as the feature
        type, location, qualifiers (e.g., product, locus_tag, external_id),

        Returns:
          a dictionary
        """
        note = {}
        if "note" in seq_feature.qualifiers:
            temp_note = ";".join(seq_feature.qualifiers["note"])
            note = note | self._process_feature_note_protein_exceptions(temp_note)
            note = note | self._process_go(temp_note)
        return note

    def _add_feature(
        self,
        seq_record,
        seq_feature,
        assembly_id,
        locus_id,
        keep_sequence,
    ):
        """
        Adds a sequence feature to the child SocialGene object

        Args:
          seq_record: The `seq_record` parameter is a BioPython `SeqRecord` object, which represents a
        biological sequence record. It contains information about the sequence, such as its sequence
        data, annotations, and features.
          seq_feature: The `seq_feature` parameter is a SeqFeature object that represents a feature
        (e.g., gene, CDS, protein) in a sequence record. It contains information such as the feature
        type, location, qualifiers (e.g., product, locus_tag, external_id)
          assembly_id: The `assembly_id` parameter is used to identify the assembly to which the
        sequence feature belongs.
          locus_id: The `locus_id` parameter is used to identify a specific locus.
          keep_sequence: The `keep_sequence` parameter is a boolean flag that determines whether the
        sequence of the added feature should be retained or not.
        """
        # Grab the locus/protein description if it exists
        if "product" in seq_feature.qualifiers:
            description = ";".join(seq_feature.qualifiers["product"]).strip()
        else:
            description = None
        if seq_feature.type == "source":
            # TODO: should this overwrite (current) or compare every iteration per assembly?
            self.assemblies[assembly_id].metadata.update(seq_feature.qualifiers)
            self.assemblies[assembly_id].loci[locus_id].metadata.update(
                seq_feature.qualifiers
            )
        if seq_feature.type in ["protein", "CDS"]:
            try:
                locus_tag = None
                external_id = None
                if "locus_tag" in seq_feature.qualifiers:
                    locus_tag = ";".join(seq_feature.qualifiers["locus_tag"])

                if "protein_id" in seq_feature.qualifiers:
                    external_id = ";".join(seq_feature.qualifiers["protein_id"])
                elif locus_tag:
                    # use locus tag as protein instead
                    external_id = locus_tag
                else:
                    # else assign random id
                    external_id = str(uuid4())
                if "translation" in seq_feature.qualifiers:
                    translation = seq_feature.qualifiers["translation"][0]
                elif "pseudo" or "pseudogene" in seq_feature.qualifiers:
                    # removes biopython potential deprecation warning to prepend Ns if non-modulo 3 seqs,
                    # may be brittle if biopython changes
                    translation = str(
                        Seq.Seq(seq_feature.extract(seq_record).seq).translate()
                    )
                    description = f"pseudo_{description}"
                else:
                    raise ValueError(
                        f"Panic!!! Not a protein or pseudo protein: {seq_feature.qualifiers}"
                    )
                uid = self.add_protein(
                    description=description,
                    external_id=external_id.strip(),
                    sequence=translation.strip(),
                    return_uid=True,
                )
                if not keep_sequence:
                    self.proteins[uid].sequence = None
                self.assemblies[assembly_id].loci[locus_id].add_feature(
                    type=seq_feature.type,
                    uid=uid,
                    external_id=external_id.strip(),
                    start=get_seqio_start(seq_feature),
                    end=get_seqio_end(seq_feature),
                    strand=seq_feature.location.strand,
                    locus_tag=locus_tag,
                    description=description,
                    **self._process_feature_note(seq_feature),
                )

            except Exception as e:
                raise e

    def _parse_record(self, seq_record, assembly_id, keep_sequence):
        """
        Parses a sequence record, adds a locus to the record, adds features
        to the record, and updates the count of loci in the file.

        Args:
          seq_record: The `seq_record` parameter is an instance of the `SeqRecord` class, which
        represents a biological sequence record. It contains information about the sequence, such as its
        sequence data, annotations, and features.
          assembly_id: The `assembly_id` parameter is an identifier for the assembly to which the
        sequence record belongs. It is used to associate the sequence record and its features with a
        specific assembly.
          keep_sequence: The parameter "keep_sequence" is a boolean value that determines whether or not
        to keep the sequence information when adding a feature.
        """
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
        **kwargs,
    ):
        """Internal function that parses a genbank file, adding assembly and taxon information,
         and then calls another function to parse the record and its features.
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
        """
        The `_open_genbank` function is used to parse a genbank file, which can handle .gz and tar
        archives.

        Args:
          input_path: The `input_path` parameter is a string that represents the path to the genbank
        file that you want to parse. It should be able to handle both uncompressed files and files
        compressed with gzip or tar. (tar hasn't been tested in a while)
        """
        # self._count_loci_in_file is mostly for sending a summary of parsed features to the logger
        log.info(f"Parsing: {input_path}")
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
            with fh.open_read(input_path) as handle:
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
