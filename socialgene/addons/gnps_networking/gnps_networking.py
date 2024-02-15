import argparse
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd

from socialgene.addons.gnps_library.gnps_library import GnpsLibrarySpectrum, GnpsLibrarySpectrumNode
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.neo4j.neo4j_element import Node, Relationship
from socialgene.utils.logging import log


def extract_parameter(xml_string):
    root = ET.fromstring(xml_string)
    return {root.attrib["name"]: root.text}


def capture_assembly_id(s, regex):
    match regex:
        case "_":
            pattern = r"^([__]+)"
        case "__":
            pattern = r"^([^_]+)__"
        case _:
            pattern = regex
    match = re.search(pattern, s)
    if match:
        return match.group(1)
    else:
        return s


# Node classes


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


class GNPS_SNETS:
    """Parses GNPS molecular network results and adds them to a Neo4j database"""

    def __init__(
        self,
        gnps_dirpath,
    ):
        self.gnps_dirpath = Path(gnps_dirpath)
        self.params_xml_path = None
        self.specnets_path = None
        self.selfloop_path = None
        self.clustersummary_path = None
        self.clusterinfo_path = None
        self.params = []
        self.filemap = {}
        self.clusterinfo_df = pd.DataFrame()
        self.specnets_df = pd.DataFrame()
        self.selfloop_df = pd.DataFrame()
        self.clustersummary_df = pd.DataFrame()
        self.workflow_uuid = str(uuid4())
        self._get_gnps_paths()
        self._get_gnps_params()
        self._file_mapping()
        self._workflow_uuid()
        self._parse_clusterinfo()
        self._parse_dfs()

    def _get_gnps_paths(self):
        """Parses the downloaded GNPS results for the necessary files, must be unzipped"""
        params_xml_path = Path(self.gnps_dirpath, "params.xml")
        specnets_path = list(self.gnps_dirpath.glob("**/result_specnets_DB/*.tsv"))
        clustersummary_path = list(
            self.gnps_dirpath.glob(
                "**/clusterinfosummarygroup_attributes_withIDs_withcomponentID/*.clustersummary"
            )
        )
        selfloop_path = list(
            self.gnps_dirpath.glob("**/networkedges_selfloop/*.selfloop")
        )
        clusterinfo_path = list(self.gnps_dirpath.glob("**/clusterinfo/*.clusterinfo"))
        for i in [specnets_path, clustersummary_path, selfloop_path, clusterinfo_path]:
            if len(i) > 1:
                raise ValueError(f"More than one file found for {i}")
            if len(i) == 0:
                raise ValueError(f"No file found for {i}")
        specnets_path = specnets_path[0]
        clustersummary_path = clustersummary_path[0]
        selfloop_path = selfloop_path[0]
        clusterinfo_path = clusterinfo_path[0]
        for i in [
            specnets_path,
            clustersummary_path,
            selfloop_path,
            params_xml_path,
            clusterinfo_path,
        ]:
            if not i.exists():
                raise ValueError(f"File not found for {i}")
        self.params_xml_path = params_xml_path
        self.specnets_path = specnets_path
        self.selfloop_path = selfloop_path
        self.clustersummary_path = clustersummary_path
        self.clusterinfo_path = clusterinfo_path

    def _get_gnps_params(self):
        """Makes a list of dicts of the GNPS run parameters from the params.xml file"""
        with open(self.params_xml_path, "r") as f:
            for i in f:
                try:
                    self.params.append(extract_parameter(i))
                except Exception:
                    continue

    def _file_mapping(self):
        """Parses the file mapping from the GNPS run parameters and creates a dict of the original filename to the mapped filename"""
        for i in self.params:
            for k, v in i.items():
                if k == "upload_file_mapping":
                    if (
                        v.split("|")[1].split("/")[0]
                        == [i["user"] for i in self.params if "user" in i][0]
                    ):
                        self.filemap[v.split("|")[0]] = v.split("|")[1].split("/")[-1]

    def _parse_clusterinfo(self):
        """Parses the clusterinfo file and replaces the mapped filename with the original filename"""
        clusterinfo_df = pd.read_csv(self.clusterinfo_path, sep="\t").replace(
            np.nan, None
        )
        clusterinfo_df.columns = [
            name.lower()
            .replace("#", "")
            .replace(" ", "_")
            .replace("(", "")
            .replace(")", "")
            .replace("/", "_")
            .replace("-", "_")
            for name in clusterinfo_df.columns
        ]
        # replaces the mapped filename with the original filename
        clusterinfo_df["original_filename"] = clusterinfo_df["filename"].apply(
            lambda x: self.filemap[x.split("/")[-1]]
        )
        self.clusterinfo_df = clusterinfo_df

    def _workflow_uuid(self):
        """Parses the workflow uuid from the params.xml file"""
        with open(self.params_xml_path, "r") as f:
            for i in f:
                if 'parameter name="uuid"' in i:
                    self.workflow_uuid = i.split(">")[1].split("<")[0]

    def _parse_dfs(self):
        """Parses the specnets, selfloop, and clustersummary files into pandas dataframes and sanitizes the column names with lowercase and underscores"""
        for i in ["specnets_path", "selfloop_path", "clustersummary_path"]:
            df = pd.read_csv(getattr(self, i), sep="\t").replace(np.nan, None)
            df.columns = [
                name.lower()
                .replace("#", "")
                .replace(" ", "_")
                .replace("(", "")
                .replace(")", "")
                .replace("/", "_")
                .replace("-", "_")
                for name in df.columns
            ]
            setattr(self, f"{i.split('_')[0]}_df", df)
        self.clustersummary_df.drop(columns=["uniquefilesources"], inplace=True)
        self.specnets_df.drop(
            columns=[
                "updateworkflowname",
                "scan",
                "spectrumfile",
                "libraryname",
                "filescanuniqueid",
                "numberhits",
                "tags",
            ],
            inplace=True,
        )

    def _filename_to_assembly(self, regex="_"):
        """Parses the assembly from the original MS filename and adds it to the clusterinfo dataframe (Assumes file is named like 'assembly_*.mzXML' or 'assembly_*.mzML' etc.)
        You will need to modify this function if your files are named differently. These assembly names must be identical to the assembly names in the Neo4j database.
        Also means all MS files must be named with a consistent pattern.

        Args:
            regex (str, optional): Pattern of MS file names. If '_' then it expects 'assembly-id_*.mzML'; if '__' then 'assembly-id__*.mzML'; can also be a custom regex pattern. Defaults to "_".
        """
        self.clusterinfo_df["assembly"] = self.clusterinfo_df[
            "original_filename"
        ].apply(lambda x: capture_assembly_id(x, regex=regex))

    def _check_db_for_assemblies(self):
        """Checks the Neo4j database for assemblies that can be linked to the GNPS results"""
        assemblies = set(self.clusterinfo_df["assembly"].unique())
        with GraphDriver() as db:
            results = db.run(
                """
                WITH $assemblies as assemblies
                UNWIND assemblies as assembly
                MATCH (a1:assembly {uid: assembly})
                RETURN assembly
                """,
                assemblies=list(assemblies),
            ).value()
        if results:
            log.warning(
                f"Assemblies from GNPS results found in db: {len(results)} of {len(assemblies)}"
            )
            log.warning(
                f"Assemblies from GNPS results not found in db: {assemblies - set(results)}"
            )

    def add_gnps_library_spectra(self):
        """Adds GNPS library hits to the Neo4j database"""
        log.info("Adding GNPS library hits to db")
        # This is potentially slow w/ a lot of trips to Neo4j but is more standardized so leaving it for now
        for i in self.specnets_df.to_dict("records"):
            a = GnpsLibrarySpectrum(**i)
            a.add_cmpd()
            # a.add_node_to_neo4j()

    def add_gnps_library_spectrum_classifications(self):
        """Adds GNPS library classifications to the Neo4j database"""
        log.info("Adding GNPS library classifications to db")
        with GraphDriver() as db:
            _ = db.run(
                """
                WITH $df as df
                UNWIND df as row
                WITH row
                WHERE row.npclassifier_class IS NOT NULL
                MATCH (gls:gnps_library_spectrum {uid: row.spectrumid})

                MERGE (npclassifier_class:npclassifier_class {uid: row.npclassifier_class})
                MERGE (npclassifier_superclass:npclassifier_superclass {uid: row.npclassifier_superclass})
                MERGE (npclassifier_pathway:npclassifier_pathway {uid: row.npclassifier_pathway})
                MERGE (npclassifier_class)-[:IS_A]->(npclassifier_superclass)
                MERGE (npclassifier_superclass)-[:IS_A]->(npclassifier_pathway)
                MERGE (gls)-[:IS_A]->(npclassifier_class)
                """,
                df=self.specnets_df.to_dict("records"),
            ).value()

    def add_gnps_library_spectrum_ion_sources(self):
        """Adds GNPS library ion sources to the Neo4j database"""
        log.info("Adding GNPS library ion sources to db")
        with GraphDriver() as db:
            _ = db.run(
                """
                WITH $df as df
                UNWIND df as row
                WITH row
                WHERE row.ion_source IS NOT NULL
                MATCH (gls:gnps_library_spectrum {uid: row.spectrumid})
                MERGE (ion_source:ion_source {uid: row.ion_source})
                MERGE (gls)-[:FROM]->(ion_source)
                """,
                df=self.specnets_df.to_dict("records"),
            ).value()

    def add_gnps_library_spectrum_instruments(self):
        """Adds GNPS library instruments to the Neo4j database"""
        log.info("Adding GNPS library instruments to db")
        with GraphDriver() as db:
            _ = db.run(
                """
                WITH $df as df
                UNWIND df as row
                WITH row
                WHERE row.instrument IS NOT NULL
                MATCH (gls:gnps_library_spectrum {uid: row.spectrumid})
                MERGE (instrument:instrument {uid: row.instrument})
                MERGE (gls)-[:FROM]->(instrument)
                """,
                df=self.specnets_df.to_dict("records"),
            ).value()

    def add_gnps_library_spectrum_organisms(self):
        """Adds GNPS library organisms to the Neo4j database"""
        log.info("Adding GNPS library organisms to db")
        with GraphDriver() as db:
            _ = db.run(
                """
                WITH $df as df
                UNWIND df as row
                WITH row
                WHERE row.organism IS NOT NULL
                MATCH (gls:gnps_library_spectrum {uid: row.spectrumid})
                MERGE (organism:organism {uid: row.organism})
                MERGE (gls)-[:FROM]->(organism)
                """,
                df=self.specnets_df.to_dict("records"),
            ).value()

    def add_gnps_clusters(self):
        """Adds GNPS networking clusters to the Neo4j database and links them to the GNPS library hits"""
        log.info("Adding GNPS clusters to db")
        with GraphDriver() as db:
            _ = db.run(
                """
                    WITH $df as df
                    UNWIND df as row
                    MERGE (c:cluster {uid: row.cluster_index, workflow_uuid: $workflow_uuid})
                        ON CREATE SET   c.defaultgroups = row.defaultgroups,
                                        c.g1 = row.g1,
                                        c.g2 = row.g2,
                                        c.g3 = row.g3,
                                        c.g4 = row.g4,
                                        c.g5 = row.g5,
                                        c.g6 = row.g6,
                                        c.gnpslinkout_cluster = row.gnpslinkout_cluster,
                                        c.gnpslinkout_network = row.gnpslinkout_network,
                                        c.mqscore = row.mqscore,
                                        c.mzerrorppm = row.mzerrorppm,
                                        c.massdiff = row.massdiff,
                                        c.rtmean = row.rtmean,
                                        c.rtmean_min = row.rtmean_min,
                                        c.rtstderr = row.rtstderr,
                                        c.uniquefilesources = row.uniquefilesources,
                                        c.uniquefilesourcescount = row.uniquefilesourcescount,
                                        c.cluster_index = row.cluster_index,
                                        c.componentindex = row.componentindex,
                                        c.number_of_spectra = row.number_of_spectra,
                                        c.parent_mass = row.parent_mass,
                                        c.precursor_charge = row.precursor_charge,
                                        c.precursor_mass = row.precursor_mass,
                                        c.sumprecursor_intensity = row.sumprecursor_intensity
                    WITH c, row
                    WHERE row.spectrumid STARTS WITH "CCMSLIB"
                    MERGE (gls:gnps_library_spectrum {uid: row.spectrumid})
                        ON CREATE SET gls.smiles = row.smiles,
                                        gls.libraryid = row.libraryid
                    WITH row, c, gls
                    MATCH (gls:gnps_library_spectrum {uid: row.spectrumid})
                    MERGE (c)-[lh:LIBRARY_HIT]->(gls)
                    SET c:gnps_library_hit
                    """,
                df=self.clustersummary_df.to_dict("records"),
                workflow_uuid=self.workflow_uuid,
            ).value()

    def add_links_between_gnps_clusters_and_assemblies(self):
        """Adds links between GNPS clusters and genome assemblies to the Neo4j database"""
        log.info("Adding links between GNPS clusters and assemblies")
        df = (
            self.clusterinfo_df.groupby(["assembly"])["clusteridx"]
            .apply(set)
            .apply(list)
        )
        for i in df.reset_index().itertuples(index=False):
            with GraphDriver() as db:
                _ = db.run(
                    """
                    WITH $row as row
                    MATCH (a1:assembly {uid: row.assembly})
                    UNWIND row.clusteridx as clusteridx
                    MATCH (c:cluster {uid: clusteridx, workflow_uuid: $workflow_uuid})
                    MERGE (c)<-[:PRODUCES]-(a1)
                    """,
                    row={"assembly": i.assembly, "clusteridx": i.clusteridx},
                    workflow_uuid=self.workflow_uuid,
                ).value()

    def add_input_spectra_to_gnps_clusters(self):
        """Adds and links input spectra to GNPS clusters in the Neo4j database"""
        log.info("Adding and linking input spectra to GNPS clusters")
        with GraphDriver() as db:
            _ = db.run(
                """
                WITH $df as df
                UNWIND df as row
                MATCH (c:cluster {uid: row.clusteridx, workflow_uuid: $workflow_uuid})
                MERGE (s:spectrum {uid: row.scan, original_filename: row.original_filename})
                    ON CREATE SET
                    s.parentmass = row.parentmass,
                    s.charge = row.charge,
                    s.rettime = row.rettime
                    s.assembly = row.assembly
                MERGE (c)-[:CONTAINS]->(s)
                """,
                df=self.clusterinfo_df.to_dict("records"),
                workflow_uuid=self.workflow_uuid,
            ).value()

    def add_links_between_gnps_clusters(self):
        """Adds links between GNPS clusters in the Neo4j database (i.e. the molecular network)"""
        log.info("Adding links between GNPS clusters")
        with GraphDriver() as db:
            _ = db.run(
                """
                WITH $df as df
                UNWIND df as row
                MATCH (c1:cluster {uid: row.clusterid1, workflow_uuid: $workflow_uuid})
                MATCH (c2:cluster {uid: row.clusterid2, workflow_uuid: $workflow_uuid})
                MERGE (c1)-[s:LINK]->(c2)
                    ON CREATE SET s.deltamz = row.deltamz,
                                  s.meh = row.meh,
                                  s.cosine = row.cosine,
                                  s.otherscore = row.otherscore,
                                  s.componentindex = row.componentindex,
                                  s.edgeannotation = row.edgeannotation
                """,
                df=self.selfloop_df.to_dict("records"),
                workflow_uuid=self.workflow_uuid,
            ).value()


def create_arg_parser():
    """ "Creates and returns the ArgumentParser object."""
    parser = argparse.ArgumentParser(
        description="Integrate a GNPS molecular network into a SocialGene Neo4j Database"
    )
    parser.add_argument(
        "--inputdir",
        help="Path to unzipped directory of GNPS molecular network download",
        default=False,
        required=True,
    )
    parser.add_argument(
        "--regex",
        help="Pattern of MS file names. If '_' then it expects 'assembly-id_*.mzML'; if '__' then 'assembly-id__*.mzML'; can also be a custom regex pattern. Defaults to '_'.",
        default="_",
    )
    return parser


def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    if not Path(args.inputdir).exists():
        raise FileNotFoundError(f"Path {args.inputdir} does not exist")
    if not Path(args.inputdir).is_dir():
        raise ValueError(f"Path {args.inputdir} is not a directory")
    gnps = GNPS_SNETS(gnps_dirpath=args.inputdir)
    gnps._filename_to_assembly(regex=args.regex)
    gnps.add_gnps_library_spectra()
    gnps.add_gnps_library_spectrum_classifications()
    gnps.add_gnps_library_spectrum_ion_sources()
    gnps.add_gnps_library_spectrum_instruments()
    gnps.add_gnps_library_spectrum_organisms()
    gnps.add_gnps_clusters()
    gnps.add_links_between_gnps_clusters_and_assemblies()
    gnps.add_input_spectra_to_gnps_clusters()
    gnps.add_links_between_gnps_clusters()


if __name__ == "__main__":
    main()
