from types import GeneratorType

from socialgene.compare_proteins.diamond import DiamondBlastp
from socialgene.compare_proteins.hmmer import CompareDomains
from socialgene.compare_proteins.mmseqs import MMseqsEasySearch


class BGCComparison:
    def __init__(self, tool):
        match tool:
            case "blastp":
                self.tool = DiamondBlastp()
            case "mmseqs2":
                self.tool = MMseqsEasySearch()
            case "hmmer":
                self.tool = CompareDomains()
            case _:
                raise ValueError(
                    f"Tool {tool} not recognised, must be one of 'blastp', 'mmseqs2', 'hmmer'"
                )

    def compare(self, bgc1, bgc2, **kwargs):
        if isinstance(bgc1, GeneratorType):
            bgc1 = list(bgc1)
        if isinstance(bgc2, GeneratorType):
            bgc2 = list(bgc2)
        return self.tool.reciprocal_hits(
            self.tool.compare_proteins(bgc1, bgc2, **kwargs),
            self.tool.compare_proteins(bgc2, bgc1, **kwargs),
        )
