# parse gff file into socialgene object
from pathlib import Path

from Bio import SeqIO

import socialgene.utils.file_handling as fh


class GFFParserMixin:
    def __init__(self):
        pass

    def _check_is_gff(self, input_path: str):
        if fh.guess_filetype(input_path) != "gff":
            raise ValueError("Input file is not a GFF file")

    def _has_fasta(self, input_path: str):
        with fh.open_read(input_path) as f:
            for line in f:
                if line.startswith("##FASTA"):
                    return True
        return False

    def _fasta_biopython(self, input_path: str):
        with fh.open_read(input_path) as file:
            for line in file:
                if line.startswith("##FASTA"):
                    break
            # read the remaining lines into fasta dictionary using biopython
            fasta_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
            if not fasta_dict:
                raise ValueError("No sequences found in FASTA section")
        return fasta_dict

    def parse_gff_file(self, input_path: str, keep_sequence: bool = True):
        # name of file without extension
        assembly_uid = Path(input_path).name
        assembly_uid = assembly_uid.split(".gff")[0]
        self._check_is_gff(input_path)
        self._has_fasta(input_path)
        nucleotide_sequence_dict = self._fasta_biopython(input_path)
        self.add_assembly(assembly_uid)
        with fh.open_read(input_path) as file:
            for line in file:
                if line.startswith("#") or line.startswith(">"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) != 9:
                    continue
                (
                    seq_id,
                    source,
                    feature_type,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    attributes,
                ) = parts
                if feature_type not in ["protein", "CDS", "pseudogene"]:
                    continue
                # get the translated sequence using biopython's translate after extracting the sequence
                # from the fasta file
                sequence = nucleotide_sequence_dict[seq_id][int(start) - 1 : int(end)]
                if strand == "-":
                    sequence = sequence.reverse_complement()
                sequence = sequence.translate()
                uid = self.add_protein(
                    description=None,
                    external_id=None,
                    sequence=str(sequence.seq),
                    return_uid=True,
                )
                if not keep_sequence:
                    self.proteins[uid].sequence = None
                if strand == "-":
                    strand = -1
                elif strand == "+":
                    strand = 1
                else:
                    strand = 0
                self.assemblies[assembly_uid].add_locus(external_id=seq_id)
                self.assemblies[assembly_uid].loci[seq_id].add_feature(
                    type=feature_type,
                    description=None,
                    external_id=None,
                    uid=uid,
                    start=int(start),
                    end=int(end),
                    strand=strand,
                )
