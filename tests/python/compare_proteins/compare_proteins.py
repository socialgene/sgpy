import os
from pathlib import Path

import pandas as pd

from socialgene.base.socialgene import SocialGene
from socialgene.compare_proteins.base import BlastTab_COLUMNS
from socialgene.compare_proteins.diamond import DiamondBlastp
from socialgene.compare_proteins.hmmer import CompareDomains
from socialgene.compare_proteins.mmseqs import MMseqsEasySearch

FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data", "compare_proteins")

sg_1848 = os.path.join(FIXTURE_DIR, "BGC0001848.pickle")
sg_1850 = os.path.join(FIXTURE_DIR, "BGC0001850.pickle")


def test_DiamondBlastp_compare_proteins_dataframe():
    a1 = SocialGene().eat_pickle(sg_1848)
    a2 = SocialGene().eat_pickle(sg_1850)
    z1 = DiamondBlastp()
    z1 = z1.compare_proteins(a1.proteins.values(), a2.proteins.values())
    z2 = pd.read_csv(
        Path(FIXTURE_DIR, "test_DiamondBlastp.csv"), dtype=BlastTab_COLUMNS
    )
    pd.testing.assert_frame_equal(z1, z2, check_names=False)


def test_MMseqsEasySearch_compare_proteins_dataframe():
    a1 = SocialGene().eat_pickle(sg_1848)
    a2 = SocialGene().eat_pickle(sg_1850)
    z1 = MMseqsEasySearch()
    z1 = z1.compare_proteins(a1.proteins.values(), a2.proteins.values())
    z2 = pd.read_csv(
        Path(FIXTURE_DIR, "test_MMseqsEasySearch.csv"), dtype=BlastTab_COLUMNS
    )
    pd.testing.assert_frame_equal(z1, z2, check_names=False)


def test_CompareDomains_compare_proteins_dataframe():
    a1 = SocialGene().eat_pickle(sg_1848)
    a2 = SocialGene().eat_pickle(sg_1850)
    z1 = CompareDomains()
    z1 = z1.compare_proteins(a1.proteins.values(), a2.proteins.values())
    z2 = pd.read_csv(Path(FIXTURE_DIR, "test_CompareDomains.csv"))
    pd.testing.assert_frame_equal(z1, z2, check_names=False)
