# python dependencies
import base64
import hashlib
import argparse

# external dependencies

# internal dependencies


parser = argparse.ArgumentParser(description="Hash a string/amino acids")
parser.add_argument(
    "--input",
    metavar="string",
    help="single string to hash",
    required=True,
)


def sha512t24u_hasher(input):
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7714221/
    # To standardize to caps-only input, use hash_aminos()
    data = bytes(input, "utf8")
    sha512_digest = hashlib.sha512(data).digest()[:24]
    sha512t24u = base64.urlsafe_b64encode(sha512_digest).decode("ascii")
    return sha512t24u


def hash_aminos(input):
    # make sure everything is uppercase before hashing
    return sha512t24u_hasher(input.upper())
