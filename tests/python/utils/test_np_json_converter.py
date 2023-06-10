import datetime

import numpy as np

from socialgene.utils.np_json_converter import np_json_converter

test_data = np.array([0], dtype=int)


def test_np_json_converter_int():
    assert isinstance(np_json_converter(np.int_(1)), int)


def test_np_json_converter_float():
    assert isinstance(np_json_converter(np.double(1)), float)


def test_np_json_converter_array():
    assert isinstance(np_json_converter(np.array([0], dtype=int)), list)


def test_np_json_converter_datetime():
    assert isinstance(np_json_converter(datetime.datetime.now()), str)
