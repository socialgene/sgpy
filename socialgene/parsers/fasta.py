from functools import partial
from io import StringIO
from pathlib import Path
from uuid import uuid4

from Bio import SeqIO

import socialgene.utils.file_handling as fh
from socialgene.utils.logging import log

fasta_prefixes = {
    "swissprot": "sp|",
    "trembl": "tr|",
    # "uniclust1": "uc",
    # "uniclust2": "cl|",
    # "genbank": "gb|",
    # "ncbi_ref": "ref|",
    # "patents": "pat|",
    # "ncbi_gi": "gi|",
}


def get_after_pipe(x, n=1):
    return x.split("|")[n]


fasta_prefixes_get_identifier = {
    "swissprot": partial(get_after_pipe, n=1),
    "trembl": partial(get_after_pipe, n=1),
}


def check_prefix(id):
    for k, v in fasta_prefixes.items():
        if id.startswith(v):
            return k
    return None


class FastaParserMixin:
    def __init__(self):
        pass

    @staticmethod
    def parse_id(id):
        try:
            prefix_key = check_prefix(id)
            if prefix_key:
                return fasta_prefixes_get_identifier[prefix_key](id)
        except Exception:
            return id
        return id

    # havent tested/worked on this since lots of updates
    def parse_fasta_file(
        self,
        input: str = None,
        defline_magic: bool = False,
    ):
        """Parse a protein fasta file

        Args:
            input (str): path to fasta file
            defline_magic (bool): Parse out fasta deflines with defined structures, currently only works for uniprot (deflines that begin with: sp| or tr|)
        """
        # if not appending, reset self.records
        if isinstance(input, str) and input.startswith(">"):
            assembly_id = str(uuid4())
            _open = StringIO(input)
            input_from = "input string"
        else:
            input = Path(input)
            try:
                assembly_id = input.stem
            except Exception:
                assembly_id = str(uuid4())
            _open = fh.open_read(input)
            input_from = input
        with _open as handle:
            self.add_assembly(uid=assembly_id, parent=self)
            self.assemblies[assembly_id].add_locus(external_id=assembly_id)
            record_counter = 0
            for seq_record in SeqIO.parse(handle, "fasta"):
                # probably could use some more defensive programming here
                if defline_magic:
                    id = self.parse_id(seq_record.id)
                else:
                    id = seq_record.id
                uid = self.add_protein(
                    description=seq_record.description,
                    external_id=id,
                    sequence=str(seq_record.seq),
                    return_uid=True,
                )
                self.assemblies[assembly_id].loci[assembly_id].add_feature(
                    type="protein",
                    description=seq_record.description,
                    external_id=id,
                    uid=uid,
                    start=record_counter,
                    end=record_counter,
                    strand=0,
                )
                record_counter += 1
        log.info(f"Read {record_counter} proteins from {input_from}")
