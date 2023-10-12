from collections.abc import Generator
from typing import List

import numpy as np

from socialgene.utils.logging import log


def chunk_a_list_with_numpy(input_list: List, n_chunks: int) -> Generator:
    """Chunk a list into n-lists

    Args:
        input_list (list): input list to be chunked
        n_chunks (int): integer of number chunks to split list into
    """
    if len(input_list) < n_chunks:
        log.info(
            "chunk_a_list_with_numpy(): n_chunks is < len(input_list), proceeding with len(input_list)"
        )
        n_chunks = len(input_list)
    return (list(i) for i in np.array_split(input_list, n_chunks))
