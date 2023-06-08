from socialgene.utils.ping import ping, get_size
import urllib
import pytest


def test_ping_ftp_nih():
    assert ping("https://ftp.ncbi.nlm.nih.gov/robots.txt") == 200


def test_ping_fail():
    with pytest.raises(urllib.error.HTTPError):
        ping("https://ftp.ncbi.nlm.nih.gov/robots.txt2")


def test_ping_1():
    assert get_size("https://ftp.ncbi.nlm.nih.gov/robots.txt", 26) == 26


def test_ping_2():
    assert get_size("https://ftp.ncbi.nlm.nih.gov/robots.txt", 1) == 26


def test_ping_3():
    with pytest.raises(ValueError):
        get_size("https://ftp.ncbi.nlm.nih.gov/robots.txt", 2000)
