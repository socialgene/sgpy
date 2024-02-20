import re
import xml.etree.ElementTree as ET
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd
from socialgene.addons.chemistry.nr import ChemicalCompoundNode
from socialgene.addons.gnps_library.nr import GnpsLibrarySpectrumNode, GnpsLibrarySpectrumToNPClassifierClass, GnpsLibrarySpectrumToNPClassifierPathway, GnpsLibrarySpectrumToNPClassifierSuperclass, GnpsLibraryToChem
from socialgene.addons.gnps_library.parse import GNPSLibrarySpectrum
from socialgene.addons.gnps_networking.nr import ClusterNode, ClusterToAssembly, LibraryHitRel, MolecularNetwork
from socialgene.addons.npclassifier.nr import NPClassifierClass, NPClassifierPathway, NPClassifierSuperclass
from socialgene.base.chem import ChemicalCompound

from socialgene.neo4j.neo4j import GraphDriver
from socialgene.nextflow.nodes import ASSEMBLY
from socialgene.utils.logging import log


def extract_parameter(xml_string):
    root = ET.fromstring(xml_string)
    return {root.attrib["name"]: root.text}


def capture_assembly_id(s, regex):
    match regex:
        case "_":
            pattern = r"^[^_]*"
        case "__":
            pattern = r"^([^_]+)__"
        case _:
            pattern = regex
    match = re.search(pattern, s)
    if match:
        return match.group(0)
    else:
        return s




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
        self.library_nodes=set()
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

    @property
    def library_hit_nodes(self):
        """Creates GNPS library hit nodes"""
        for i in self.specnets_df.to_dict("records"):
            try:
                a = GnpsLibrarySpectrumNode()
                a.fill_from_dict(i)
                a.properties['uid'] = str(i["spectrumid"])
                self.library_nodes.add(a)
            except Exception as e:
                log.debug(f"Error creating library hit node for {i['spectrumid']}: {e}")

    def add_library_hit_nodes_to_neo4j(self, create=False):
        """Adds GNPS library hit nodes to the Neo4j database"""
        list(self.library_nodes)[0].add_multiple_to_neo4j(list(self.library_nodes), create=create)


    def link_npclassifiers(self):
        """Classifies GNPS library hits by linking them to NPClassifier nodes"""
        to_loop = [("npclassifier_superclass", NPClassifierSuperclass, GnpsLibrarySpectrumToNPClassifierSuperclass),
        ("npclassifier_class", NPClassifierClass, GnpsLibrarySpectrumToNPClassifierClass),
        ("npclassifier_pathway", NPClassifierPathway, GnpsLibrarySpectrumToNPClassifierPathway)]
        for name,nodeclass,relclass in to_loop:
            df = self.specnets_df[~self.specnets_df[name].isnull()]
            df = df[['spectrumid', name]]
            npc_set=set()
            gnp_set=set()
            rel_set=set()
            for i in df.to_dict("records"):
                npc = nodeclass(properties={'uid': i[name]})
                gnp = GnpsLibrarySpectrumNode(properties={'uid': i['spectrumid']})
                rel = relclass(start=gnp, end=npc)
                # accumulate unique nodes and relationships to bulk add
                npc_set.add(npc)
                gnp_set.add(gnp)
                rel_set.add(rel)
            npc_set = list(npc_set)
            gnp_set = list(gnp_set)
            rel_set = list(rel_set)
            for i in [npc_set, gnp_set, rel_set]:
                i[0].add_multiple_to_neo4j(i, create=False)

    def create_cluster_nodes(self, create=False):
        clusternodes=set()
        for i in self.clustersummary_df.to_dict("records"):
            a=ClusterNode()
            d=i | {"workflow_uuid": self.workflow_uuid}
            a.fill_from_dict(d)
            clusternodes.add(a)
        clusternodes = list(clusternodes)
        clusternodes[0].add_multiple_to_neo4j(clusternodes, create=create)

    def link_cluster_to_library(self):
        df = self.clustersummary_df[['spectrumid','cluster_index']]
        df = df[~df['spectrumid'].isnull()]
        s_set=set()
        e_set=set()
        r_set=set()
        for i in df.to_dict("records"):
            s = ClusterNode(properties={'cluster_index': i['cluster_index'], 'workflow_uuid': self.workflow_uuid})
            e = GnpsLibrarySpectrumNode(properties={'uid': i['spectrumid']})
            r = LibraryHitRel(start=s, end=e)
            s_set.add(s)
            e_set.add(e)
            r_set.add(r)
        s_set = list(s_set)
        e_set = list(e_set)
        r_set = list(r_set)
        for i in [s_set, e_set, r_set]:
            i[0].add_multiple_to_neo4j(i, create=False)


    def link_network(self, create=False):
        s_set=set()
        e_set=set()
        r_set=set()
        for i in self.selfloop_df.to_dict("records"):
            s = ClusterNode(properties={'cluster_index': i['clusterid1'], 'workflow_uuid': self.workflow_uuid})
            e = ClusterNode(properties={'cluster_index': i['clusterid2'], 'workflow_uuid': self.workflow_uuid})
            r = MolecularNetwork(start=s, end=e)
            if i['edgeannotation'] == " ":
                i['edgeannotation'] = None
            r.fill_from_dict(i)
            s_set.add(s)
            e_set.add(e)
            r_set.add(r)
        s_set = list(s_set)
        e_set = list(e_set)
        r_set = list(r_set)
        for i in [s_set, e_set]:
            # create=False because nodes should already exist so we want to merge not create
            i[0].add_multiple_to_neo4j(i, create=False)
        r_set[0].add_multiple_to_neo4j(r_set, create=create)


    def link_assembly_to_cluster(self, create=False):
        df = self.clusterinfo_df[['assembly','clusteridx']]
        df = df[~df['assembly'].isnull()]
        s_set=set()
        e_set=set()
        r_set=set()
        for i in df.to_dict("records"):
            s = ClusterNode(properties={'cluster_index': i['clusteridx'], 'workflow_uuid': self.workflow_uuid})
            e = ASSEMBLY(properties={'uid': i['assembly']})
            r = ClusterToAssembly(start=s, end=e)
            s_set.add(s)
            e_set.add(e)
            r_set.add(r)
        s_set = list(s_set)
        e_set = list(e_set)
        r_set = list(r_set)
        for i in [s_set, e_set]:
            # create=False because nodes should already exist so we want to merge not create
            i[0].add_multiple_to_neo4j(i, create=False)
        r_set[0].add_multiple_to_neo4j(r_set, create=create)


    def link_library_to_chem(self):
        all_chem_nodes = set()
        all_rels = set()
        for i in self.library_nodes:
            chemstring=None
            if "inchi" in i.properties and i.properties['inchi']:
                chemstring = i.properties['inchi']
            elif "smiles" in i.properties and i.properties['smiles']:
                chemstring = i.properties['smiles']
            if chemstring:
                try:
                    cmpd=ChemicalCompound(chemstring)
                    a=ChemicalCompoundNode()
                    a.fill_from_dict(cmpd.base_properties | cmpd.hash_dict)
                    all_chem_nodes.add(a)
                    all_rels.add(GnpsLibraryToChem(start=i, end=a))
                except Exception as e:
                    log.debug(f"Error creating chemical compound node for {i.properties['uid']}: {e}")
        all_chem_nodes = list(all_chem_nodes)
        all_rels = list(all_rels)
        all_chem_nodes[0].add_multiple_to_neo4j(all_chem_nodes, create=False)
        all_rels[0].add_multiple_to_neo4j(all_rels, create=False)


