import gzip
from functools import partial
from io import StringIO
from mimetypes import guess_type
from pathlib import Path
from uuid import uuid4

from Bio import SeqIO

from socialgene.utils.logging import log


class FastaParserMixin:
    def __init__(self):
        pass

    # havent tested/worked on this since lots of updates
    def parse_fasta_file(self, input=None):
        """Parse a protein fasta file

        Args:
            input (str): path to fasta file
        """
        # if not appending, reset self.records
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
            self.add_assembly(uid=assembly_id)
            self.assemblies[assembly_id].add_locus(external_id=assembly_id)
            record_counter = 0
            for seq_record in SeqIO.parse(handle, "fasta"):
                # probably could use some more defensive programming here
                hash_id = self.add_protein(
                    description=seq_record.description,
                    external_protein_id=seq_record.id,
                    sequence=str(seq_record.seq),
                    return_uid=True,
                )
                # set the start/end to the order of records
                self.assemblies[assembly_id].loci[assembly_id].add_feature(
                    type="protein",
                    description=seq_record.description,
                    protein_id=seq_record.id,
                    protein_hash=hash_id,
                    start=record_counter,
                    end=record_counter,
                    strand=0,
                )
                record_counter += 1
        log.info(f"Read {record_counter} proteins from {input}")

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
            self.add_assembly(uid=assembly_id)
            self.assemblies[assembly_id].add_locus(external_id=assembly_id)
            for seq_record in SeqIO.parse(handle, "fasta"):
                prot_hash = self.add_protein(
                    description=seq_record.description,
                    external_protein_id=seq_record.id,
                    sequence=str(seq_record.seq),
                    return_uid=True,
                )
                self.assemblies[assembly_id].loci[assembly_id].add_feature(
                    type="CDS",
                    protein_hash=prot_hash,
                )
                record_counter += 1
                count_proteins_in_file += 1
        log.info(f"Read {count_proteins_in_file} proteins from input string")
