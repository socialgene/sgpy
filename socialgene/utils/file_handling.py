import bz2
import gzip
import lzma
import shutil
import tarfile
from contextlib import contextmanager
from enum import Enum, auto
from pathlib import Path
from typing import IO

from socialgene.utils.logging import log


class Compression(Enum):
    bzip2 = auto()
    gzip = auto()
    xz = auto()
    uncompressed = auto()


def is_compressed(filepath: Path) -> Compression:
    """
    Determines the compression type of a file based on its signature.

    Args:
      filepath (Path): The `filepath` parameter is a `Path` object that represents the path to the file
    that you want to check for compression.

    Returns:
      The function `is_compressed` returns the type of compression used for the file specified by the
    `filepath` parameter. The possible return values are `Compression.gzip`, `Compression.bzip2`,
    `Compression.xz`, or `Compression.uncompressed`.
    """
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
    """
    The function `gunzip` decompresses a gzip file and saves the decompressed file with the same name
    but without the ".gz" extension.

    Args:
      filepath (Path): Path to the compressed file that you want to
    decompress.
    """
    if is_compressed(filepath).name == "gzip":
        # remove ".gz" for the new filepath
        new_path = Path(filepath.parents[0], filepath.stem)
        log.info(f"Started decompressing: {str(Path(filepath))}")
        # open gz, decompress, write back out to new file
        with gzip.open(filepath, "rb") as f_in:
            with open(new_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        filepath.unlink()
        log.info(f"Finished decompressing: {str(Path(filepath))}")
        return new_path


@contextmanager
def open_read(filepath: Path) -> IO:
    """
    The function `open_read` opens a file for reading, taking into account different compression
    formats.

    Args:
      filepath (Path): The `filepath` parameter is the path to the file that you want to open and read.
    It should be a string representing the file path.
    """
    filepath_compression = is_compressed(filepath)
    if filepath_compression == Compression.gzip:
        f = gzip.open(filepath, "rt")
    elif filepath_compression == Compression.bzip2:
        f = bz2.open(filepath, "rt")
    elif filepath_compression == Compression.xz:
        f = lzma.open(filepath, "rt")
    else:
        f = open(filepath, "rt")
    try:
        yield f
    finally:
        f.close()


@contextmanager
def open_write(filepath: str, mode="w", compression: str = None) -> IO:
    """Open a file for writing

    Args:
        filepath (str): input filepath
        mode (str, optional): modes to open file for writing (only use "a", "w", or "r"- "b" will auto-applied). Defaults to "w".
        compression (str, optional): which compression method to use ("gzip", "bzip", "xz", or None). Defaults to None.
    Yields:
        Iterator[IO]: context manager
    """
    filepath = Path(filepath)
    if mode not in ["w", "a", "r"]:
        raise ValueError
    match compression:
        case "gzip":
            _open = gzip.open(filepath.with_suffix(".gz"), f"{mode}")
        case "bzip":
            _open = bz2.open(filepath.with_suffix(".bz2"), f"{mode}")
        case "xz":
            _open = lzma.open(filepath.with_suffix(".xz"), f"{mode}")
        case None:
            _open = open(filepath, mode)
        case _:
            raise ValueError(
                f'compression variable must be "gzip", "bzip", "xz", or None; was {compression}'
            )
    try:
        yield _open
    finally:
        _open.close()


def check_if_tar(filepath):
    return tarfile.is_tarfile(filepath)


def guess_filetype(filepath):
    """Guess what type of file it is
    Args:
        filepath: file path of file to guess
    Returns:
      a string indicating the type of file. The possible return values are "genbank", "fasta", "gff", or
    "domtblout".
    """
    with open_read(filepath) as f:
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
