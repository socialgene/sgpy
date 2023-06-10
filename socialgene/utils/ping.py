import urllib.request
from contextlib import closing
from urllib.request import urlopen


def ping(url):
    return urllib.request.urlopen(url).getcode()


def get_size(url, expected_size):
    with closing(urlopen(url)) as response:
        cl = int(response.headers["Content-length"])
        if abs(expected_size - cl) > 1000:
            raise ValueError(
                f"Expected_size was {expected_size}. The response header claims it is {cl}. Believe what you will."
            )
    return cl
