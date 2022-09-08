# python dependencies
import base64
import hashlib
import argparse
import re

# external dependencies

# internal dependencies


parser = argparse.ArgumentParser(description="Hash a string/amino acids")
parser.add_argument(
    "--input",
    metavar="string",
    help="single string to hash",
    required=True,
)

mmseqs_will_remove_these = [
    r"^uc",
    r"^cl\|",
    r"^sp\|",
    r"^tr\|",
    r"^gb\|",
    r"^ref\|",
    r"^pdb\|",
    r"^bbs\|",
    r"^lcl\|",
    r"^pir\|\|",
    r"^prf\|\|",
    r"^gnl\|",
    r"^pat\|",
    r"^gi\|",
]

mmseqs_will_remove_these_re = re.compile("|".join(mmseqs_will_remove_these))


def sha512_hash(input):
    # To standardize to caps-only input, use hash_aminos()
    data = bytes(input, "utf8")
    sha512_digest = hashlib.sha512(data).digest()[:24]
    sha512t24u = base64.urlsafe_b64encode(sha512_digest).decode("ascii")
    # https://github.com/soedinglab/MMseqs2/issues/557
    # TODO: add prefix to avoid mmseqs issue, (this should only be a temp fix)
    if mmseqs_will_remove_these_re.match(sha512t24u):
        sha512t24u = f"mm{sha512t24u}"
    return sha512t24u


def hash_aminos(input):
    # make sure everything is uppercase before hashing
    return sha512_hash(input.upper())
