import argparse
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd

from socialgene.addons.gnps_library.nr import GnpsLibrarySpectrum, GnpsLibrarySpectrumNode
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.utils.logging import log


class ClusterNode(Node):
    neo4j_label = "gnps_cluster"
    description = "Represents a GNPS molecular networking cluster"
    required_properties = ["uid", "workflow_uuid"]
    property_specification = {
        "uid": str,
        "workflow_uuid": str,
        "defaultgroups": str,
        "g1": str,
        "g2": str,
        "g3": str,
        "g4": str,
        "g5": str,
        "g6": str,
        "gnpslinkout_cluster": str,
        "gnpslinkout_network": str,
        "mqscore": float,
        "mzerrorppm": float,
        "massdiff": float,
        "rtmean": float,
        "rtmean_min": float,
        "rtstderr": float,
        "uniquefilesources": str,
        "uniquefilesourcescount": int,
        "cluster_index": int,
        "componentindex": int,
        "number_of_spectra": int,
        "parent_mass": float,
        "precursor_charge": int,
        "precursor_mass": float,
        "sumprecursor_intensity": float,
    }


class SpectrumNode(Node):
    neo4j_label = "spectrum"
    description = "Represents a GNPS molecular networking spectrum"
    property_specification = {
        "uid": str,
        "original_filename": str,
        "parentmass": float,
        "charge": int,
        "rettime": float,
        "assembly": str,
    }


class LibraryHitRel(Relationship):
    neo4j_label = "LIBRARY_HIT"
    description = "Connects a GNPS cluster to a GNPS library hit"
    start_class = ClusterNode
    end_class = GnpsLibrarySpectrumNode

