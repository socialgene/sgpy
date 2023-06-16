from typing import TextIO

import bz2
import gzip
import lzma
import shutil
import tarfile
from contextlib import contextmanager
from enum import Enum, auto
from pathlib import Path

from socialgene.utils.logging import log


class Compression(Enum):
    bzip2 = auto()
    gzip = auto()
    xz = auto()
    uncompressed = auto()


def is_compressed(filepath: Path) -> Compression:
    with open(filepath, "rb") as f:
        signature = f.peek(8)[:8]
        if tuple(signature[:2]) == (0x1F, 0x8B):
            return Compression.gzip
        elif tuple(signature[:3]) == (0x42, 0x5A, 0x68):
            return Compression.bzip2
        elif tuple(signature[:7]) == (0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00, 0x00):
            return Compression.xz
        else:
            return Compression.uncompressed


def gunzip(filepath: Path) -> None:
    if is_compressed(filepath).name == "gzip":
        # remove ".gz" for the new filepath
        new_hmm_path = Path(filepath.parents[0], filepath.stem)
        log.info(f"Start decompressing: {str(Path(filepath).stem)}")
        # open gz, decompress, write back out to new file
        with gzip.open(filepath, "rb") as f_in:
            with open(new_hmm_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        log.info(f"Finish decompressing: {str(Path(filepath).stem)}")


@contextmanager
def open_file(filepath: Path) -> TextIO:
    filepath_compression = is_compressed(filepath)
    if filepath_compression == Compression.gzip:
        f = gzip.open(filepath, "rt")
    elif filepath_compression == Compression.bzip2:
        f = bz2.open(filepath, "rt")
    elif filepath_compression == Compression.xz:
        f = lzma.open(filepath, "rt")
    else:
        f = open(filepath, "r")
    try:
        yield f
    finally:
        f.close()


def check_if_tar(filepath):
    return tarfile.is_tarfile(filepath)


def guess_filetype(filepath):
    """Guess what type of file it is
    Args:
        filepath: file path of file to guess
    Returns:
        bool: true/false
    """
    with open_file(filepath) as f:
        l1 = f.readline()

    if l1.startswith("LOCUS "):
        return "genbank"
    if l1.startswith(">"):
        return "fasta"
    if l1.startswith("##gff-version"):
        return "gff"
    if (
        l1.replace(" ", "")
        == "#---fullsequence-----------------thisdomain-------------hmmcoordalicoordenvcoord\n"
    ):
        return "domtblout"
