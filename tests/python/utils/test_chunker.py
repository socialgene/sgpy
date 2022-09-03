from collections.abc import Generator
import pytest
import string
from socialgene.utils.chunker import chunk_a_list_with_numpy

one_to_ten = list(range(1, 10, 1))
alphabet = list(string.ascii_lowercase)[0:10]


class ParamClass(object):
    instances = [-1, -0.1, 0, "a"] + list(range(0, 12))


@pytest.fixture(params=ParamClass().instances)
def test_input(request):
    return request.param


def test_chunk_a_list_with_numpy_type(test_input):
    # send all failing types to "else:"
    if isinstance(test_input, int) and test_input > 0 and test_input < 10:
        tmp = chunk_a_list_with_numpy(input_list=one_to_ten, n_chunks=test_input)
        assert isinstance(tmp, Generator)
        tmp = [i for i in tmp]
        assert len(tmp) == test_input
    elif isinstance(test_input, int) and test_input > 9:
        # max return length should be length of input
        tmp = chunk_a_list_with_numpy(input_list=one_to_ten, n_chunks=test_input)
        assert isinstance(tmp, Generator)
        tmp = [i for i in tmp]
        assert len(tmp) == 9
    else:
        with pytest.raises(Exception) as e_info:
            chunk_a_list_with_numpy(input_list=one_to_ten, n_chunks=test_input)
            _ = e_info
