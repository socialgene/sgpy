import socialgene.utils.file_handling as fh
from socialgene.parsers.fasta import FastaParserMixin
from socialgene.parsers.genbank import GenbankParserMixin


class SequenceParser(GenbankParserMixin, FastaParserMixin):
    def __init__(self):
        super().__init__()

    def parse(self, filepath, **kwargs):
        """Parse sequence files (main function)

        Args:
            filepath (str): path to sequence file
        """
        # TODO: fh.check_if_tar(filepath=filepath)
        if isinstance(filepath, str) and filepath.startswith(">"):
            self.parse_fasta_file(input=filepath, **kwargs)
            return
        else:
            filetype = fh.guess_filetype(filepath)
            if filetype == "genbank":
                self._open_genbank(input_path=filepath, **kwargs)
                return
            elif filetype == "fasta":
                self.parse_fasta_file(input=filepath, **kwargs)
                return

        raise NotImplementedError(
            "May not be implemented, or you need to use the genbank/fasta parser directly. (e.g. for tar archives)"
        )
