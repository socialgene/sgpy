from collections import OrderedDict
import csv
from socialgene.config import env_vars
from socialgene.utils.logging import log
import platform
from importlib.metadata import version


class Parameters:
    __slots__ = [
        "HMMSEARCH_IEVALUE",
        "HMMSEARCH_BACKGROUND",
        "HMMSEARCH_BIASFILTER",
        "HMMSEARCH_NULL2",
        "HMMSEARCH_SEED",
        "HMMSEARCH_Z",
        "HMMSEARCH_DOMZ",
        "HMMSEARCH_F1",
        "HMMSEARCH_F2",
        "HMMSEARCH_F3",
        "HMMSEARCH_E",
        "HMMSEARCH_DOME",
        "HMMSEARCH_INCE",
        "HMMSEARCH_INCDOME",
        "HMMSEARCH_BITCUTOFFS",
        "platform",
        "architecture",
        "py_executable",
        "py_version",
        "sgpy_version",
        "genome_download_command",
        "mmseqs_levels",
        "mmseqs_args",
        "blast_args"
    ]
    def __init__(self, genome_download_command = None, mmseqs_levels = None, mmseqs_args = None, blast_args = None):
        for i in self.__slots__:
            setattr(self, i, env_vars.get(i, None))
        self.platform = platform.sys.platform
        self.architecture = " ".join(platform.architecture())
        self.py_executable = platform.sys.executable
        self.py_version = platform.sys.version
        self.sgpy_version = version('socialgene')
        self.genome_download_command = genome_download_command
        self.mmseqs_levels = mmseqs_levels
        self.mmseqs_args = mmseqs_args
        self.blast_args = blast_args
    def all_attributes(self):
        return OrderedDict({i: self.__getitem__(i) for i in self.__slots__})
    def __getitem__(self, item):
        try:
            return self.__getattribute__(item)
        except Exception as e:
            log.debug(e)
            return None
    def write_tsv(self, outpath):
        with open(outpath, "a") as f:
            writer = csv.writer(
                f,
                quotechar='"', quoting=csv.QUOTE_MINIMAL,
                delimiter="\t",
            )
            writer.writerows([list(self.all_attributes().values())])
