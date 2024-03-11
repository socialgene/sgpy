import re
import xml.etree.ElementTree as ET
from pathlib import Path
from uuid import uuid4

import numpy as np
import pandas as pd

from socialgene.addons.gnps_library.nr import (
    GnpsLibrarySpectrumNode,
    GnpsLibrarySpectrumToNPClassifierClass,
    GnpsLibrarySpectrumToNPClassifierPathway,
    GnpsLibrarySpectrumToNPClassifierSuperclass,
    GnpsLibraryToChem,
)
from socialgene.addons.gnps_networking.nr import (
    ClusterNode,
    ClusterToSpectrum,
    LibraryHitRel,
    MassSpecFileNode,
    MassSpecFileToAssembly,
    MassSpecFileToSpectrum,
    MolecularNetwork,
    SpectrumNode,
)
from socialgene.addons.npclassifier.nr import (
    NPClassifierClass,
    NPClassifierPathway,
    NPClassifierSuperclass,
)
from socialgene.base.chem import ChemicalCompound
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.nextflow.nodes import ASSEMBLY
from socialgene.utils.file_handling import open_read
from socialgene.utils.logging import log


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
        gnps_dirpath: str = None,
        params_xml_path: str = None,
        specnets_path: str = None,
        selfloop_path: str = None,
        clustersummary_path: str = None,
        clusterinfo_path: str = None,
        map_path: str = None,
    ):
        self.gnps_dirpath = gnps_dirpath
        if self.gnps_dirpath:
            self.gnps_dirpath = Path(self.gnps_dirpath)
        self.params_xml_path = params_xml_path
        self.specnets_path = specnets_path
        self.selfloop_path = selfloop_path
        self.clustersummary_path = clustersummary_path
        self.clusterinfo_path = clusterinfo_path
        self.map_path = map_path
        self.gnps_file_mapping = {}
        self.clusterinfo_df = pd.DataFrame()
        self.specnets_df = pd.DataFrame()
        self.selfloop_df = pd.DataFrame()
        self.clustersummary_df = pd.DataFrame()
        self.workflow_uuid = str(uuid4()) + "_sg_assigned"
        self.task = str(uuid4()) + "_sg_assigned"
        self.library_nodes = set()
        for i in [
            "specnets_path",
            "selfloop_path",
            "clustersummary_path",
            "clusterinfo_path",
            "map_path",
        ]:
            if getattr(self, i):
                setattr(self, i, Path(getattr(self, i)))

    def _parse(self):
        self._get_gnps_filemapping()
        self._workflow_uuid()
        self._task_id()
        self._parse_clusterinfo()
        self._parse_dfs()

    def _find_gnps_paths(self):
        """Parses the downloaded GNPS results for the necessary files, must be unzipped"""
        params_xml_path = Path(self.gnps_dirpath, "params.xml")
        specnets_path = list(self.gnps_dirpath.glob("**/result_specnets_DB/*.tsv"))
        clustersummary_path = list(
            self.gnps_dirpath.glob(
                "**/clusterinfosummarygroup_attributes_withIDs_withcomponentID/*.clustersummary"
            )
        )
        selfloop_path = list(self.gnps_dirpath.glob("**/networkedges_selfloop/*"))
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

    def _get_gnps_filemapping(self):
        with open_read(self.params_xml_path) as f:
            xml_root = ET.fromstring(f.read())
        for i in xml_root.iter("parameter"):
            if i.attrib["name"] == "upload_file_mapping":
                try:
                    self.gnps_file_mapping[i.text.split("|")[0]] = i.text.split("/")[-1]
                except Exception as e:
                    log.debug(f"Error parsing file mapping: {e}")

    def _parse_clusterinfo(self):
        """Parses the clusterinfo file and replaces the mapped filename with the original filename"""
        with open_read(self.clusterinfo_path) as f:
            clusterinfo_df = pd.read_csv(f, sep="\t").replace(np.nan, None)
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
            lambda x: self.gnps_file_mapping[x.split("/")[-1]]
        )
        self.clusterinfo_df = clusterinfo_df

    def _workflow_uuid(self):
        """Parses the workflow uuid from the params.xml file"""
        with open_read(self.params_xml_path) as f:
            for i in f:
                if 'parameter name="uuid"' in i:
                    self.workflow_uuid = i.split(">")[1].split("<")[0]

    def _task_id(self):
        """Parses the task id from the params.xml file"""
        with open_read(self.params_xml_path) as f:
            for i in f:
                if 'parameter name="task"' in i:
                    self.task = i.split(">")[1].split("<")[0]

    def _parse_dfs(self):
        """Parses the specnets, selfloop, and clustersummary files into pandas dataframes and sanitizes the column names with lowercase and underscores"""
        for i in ["specnets_path", "selfloop_path", "clustersummary_path"]:
            with open_read(getattr(self, i)) as f:
                df = pd.read_csv(f, sep="\t").replace(np.nan, None)
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
        try:
            self.clustersummary_df.drop(columns=["uniquefilesources"], inplace=True)
        except KeyError:
            log.debug("No column 'uniquefilesources' in clustersummary file")
        try:
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
        except KeyError:
            log.debug(
                f"{[i for i in ['updateworkflowname', 'scan', 'spectrumfile', 'libraryname', 'filescanuniqueid', 'numberhits', 'tags'] if i not in self.specnets_df.columns]}"
            )

    @staticmethod
    def _check_db_for_assemblies(assemblies):
        """Checks the Neo4j database for assemblies that can be linked to the GNPS results"""
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
            log.info(
                f"Assemblies in GNPS results found in db: {len(results)} of {len(assemblies)}"
            )
            log.info(
                f"Assemblies in GNPS results not found in db: {set(assemblies) - set(results)}"
            )
        else:
            log.warning(
                "No assemblies in GNPS results matched to assemblies in the database"
            )
        return results

    # Add nodes to Neo4j
    def _add_input_spectra(self, create=False):
        def temp(i):
            spec = SpectrumNode()
            d = i | {"workflow_uuid": self.workflow_uuid}
            spec.fill_from_dict(d)
            return spec

        SpectrumNode.add_multiple_to_neo4j(
            [temp(i) for i in self.clusterinfo_df.to_dict("records")], create=create
        )

    def _add_input_ms_files(self, create=False):
        def temp(i):
            spec = MassSpecFileNode()
            d = {
                "workflow_uuid": self.workflow_uuid,
                "gnps_filename": i["filename"],
                "filename": i["original_filename"],
            }
            spec.fill_from_dict(d)
            return spec

        MassSpecFileNode.add_multiple_to_neo4j(
            [temp(i) for i in self.clusterinfo_df.to_dict("records")], create=create
        )

    def _create_library_hit_nodes(self):
        """Creates GNPS library hit nodes"""
        for i in self.specnets_df.to_dict("records"):
            try:
                a = GnpsLibrarySpectrumNode()
                a.fill_from_dict(i)
                a.properties["uid"] = str(i["spectrumid"])
                self.library_nodes.add(a)
            except Exception as e:
                log.debug(f"Error creating library hit node for {i['spectrumid']}: {e}")

    def _add_library_hit_nodes(self, create=False):
        """Adds GNPS library hit nodes to the Neo4j database"""
        list(self.library_nodes)[0].add_multiple_to_neo4j(
            list(self.library_nodes), create=create
        )

    def _add_cluster_nodes(self, create=False):
        clusternodes = set()
        for i in self.clustersummary_df.to_dict("records"):
            a = ClusterNode()
            d = i | {"workflow_uuid": self.workflow_uuid, "task": self.task}
            a.fill_from_dict(d)
            clusternodes.add(a)
        clusternodes = list(clusternodes)
        clusternodes[0].add_multiple_to_neo4j(clusternodes, create=create)

    # Add rels to Neo4j

    def _link_cluster_to_spectrum(self):
        rels = set()
        for i in self.clusterinfo_df.to_dict("records"):
            cn = ClusterNode()
            cn.fill_from_dict(
                {
                    "cluster_index": i["clusteridx"],
                    "workflow_uuid": self.workflow_uuid,
                    "task": self.task,
                }
            )
            sn = SpectrumNode()
            sn.fill_from_dict(i | {"workflow_uuid": self.workflow_uuid})
            rels.add(ClusterToSpectrum(start=sn, end=cn))
        rels = list(rels)
        rels[0].add_multiple_to_neo4j(rels, create=False)

    def _validate_map(self):
        # columns: assembly_id, mass_spec_filename
        if not self.map_path or not Path(self.map_path).exists():
            raise ValueError("No map file found")
        with open_read(self.map_path) as f:
            map_df = pd.read_csv(
                f,
                sep="\t",
                header=None,
            ).replace(np.nan, None)
        # should only have two columns
        if len(map_df.columns) != 2:
            raise ValueError("Map file should only have two columns")
        map_df.columns = ["assembly", "mass_spec_file"]
        # get first column as list
        assemblies = set(map_df.iloc[:, 0])
        matched_assemblies = self._check_db_for_assemblies(assemblies)
        # filter df for matched assemblies
        map_df = map_df[map_df.iloc[:, 0].isin(matched_assemblies)]
        # change second column to Path().name
        map_df.iloc[:, 1] = map_df.iloc[:, 1].apply(lambda x: Path(x).name)
        return map_df

    def _link_ms_file_to_assembly(self):
        map_df = self._validate_map()
        # merge map_df and clusterinfo_df, merge on mass_spec_file and original_filename
        df = self.clusterinfo_df.merge(
            map_df, left_on="original_filename", right_on="mass_spec_file", how="inner"
        )
        rels = set()
        for i in df.to_dict("records"):
            cn = MassSpecFileNode()
            d = {
                "workflow_uuid": self.workflow_uuid,
                "gnps_filename": i["filename"],
                "filename": i["original_filename"],
            }
            cn.fill_from_dict(d)
            sn = ASSEMBLY()
            d = {"uid": i["assembly"]}
            sn.fill_from_dict(d)
            rels.add(MassSpecFileToAssembly(start=cn, end=sn))
        rels = list(rels)
        rels[0].add_multiple_to_neo4j(rels, create=False)

    def _link_spectrum_to_ms_file(self):
        rels = set()
        for i in self.clusterinfo_df.to_dict("records"):
            cn = MassSpecFileNode()
            d = {
                "workflow_uuid": self.workflow_uuid,
                "gnps_filename": i["filename"],
                "filename": i["original_filename"],
            }
            cn.fill_from_dict(d)
            sn = SpectrumNode()
            d = i | {"workflow_uuid": self.workflow_uuid}
            sn.fill_from_dict(d)
            rels.add(MassSpecFileToSpectrum(start=cn, end=sn))
        rels = list(rels)
        rels[0].add_multiple_to_neo4j(rels, create=False)

    def _link_npclassifiers(self):
        """Classifies GNPS library hits by linking them to NPClassifier nodes"""
        to_loop = [
            (
                "npclassifier_superclass",
                NPClassifierSuperclass,
                GnpsLibrarySpectrumToNPClassifierSuperclass,
            ),
            (
                "npclassifier_class",
                NPClassifierClass,
                GnpsLibrarySpectrumToNPClassifierClass,
            ),
            (
                "npclassifier_pathway",
                NPClassifierPathway,
                GnpsLibrarySpectrumToNPClassifierPathway,
            ),
        ]
        for name, nodeclass, relclass in to_loop:
            if name not in self.specnets_df.columns:
                # not all GNPS results will have npclassifier results
                continue
            df = self.specnets_df[~self.specnets_df[name].isnull()]
            df = df[["spectrumid", name]]
            npc_set = set()
            gnp_set = set()
            rel_set = set()
            for i in df.to_dict("records"):
                npc = nodeclass(properties={"uid": i[name]})
                gnp = GnpsLibrarySpectrumNode(properties={"uid": i["spectrumid"]})
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

    def _link_cluster_to_library_no_spec_id(self):
        if "spectrumid" not in self.clustersummary_df.columns:
            if "libraryid" not in self.clustersummary_df.columns:
                log.warn("No libraryid or spectrumid in clustersummary file")
                return
        df = self.clustersummary_df[["libraryid", "cluster_index"]]
        df = df[~df["libraryid"].isnull()]
        s_set = set()
        e_set = set()
        r_set = set()
        for i in df.to_dict("records"):
            s = ClusterNode(
                properties={
                    "cluster_index": i["cluster_index"],
                    "workflow_uuid": self.workflow_uuid,
                    "task": self.task,
                }
            )
            e = GnpsLibrarySpectrumNode(properties={"uid": i["libraryid"]})
            r = LibraryHitRel(start=s, end=e)
            s_set.add(s)
            e_set.add(e)
            r_set.add(r)
        s_set = list(s_set)
        e_set = list(e_set)
        r_set = list(r_set)
        for i in [s_set, e_set, r_set]:
            i[0].add_multiple_to_neo4j(i, create=False)

    def _link_cluster_to_library(self):
        if "spectrumid" not in self.clustersummary_df.columns:
            # not all GNPS results will have library hits
            log.warn(
                "Older versions of GNPS results may have library ids but not spectrumid? Skipping."
            )
            try:
                self._link_cluster_to_library_no_spec_id()
            except Exception as e:
                raise ValueError(f"Failed to link clusters to library hits: {e}")
            return
        df = self.clustersummary_df[["spectrumid", "cluster_index"]]
        df = df[~df["spectrumid"].isnull()]
        s_set = set()
        e_set = set()
        r_set = set()
        for i in df.to_dict("records"):
            s = ClusterNode(
                properties={
                    "cluster_index": i["cluster_index"],
                    "workflow_uuid": self.workflow_uuid,
                    "task": self.task,
                }
            )
            e = GnpsLibrarySpectrumNode(properties={"uid": i["spectrumid"]})
            r = LibraryHitRel(start=s, end=e)
            s_set.add(s)
            e_set.add(e)
            r_set.add(r)
        s_set = list(s_set)
        e_set = list(e_set)
        r_set = list(r_set)
        for i in [s_set, e_set, r_set]:
            i[0].add_multiple_to_neo4j(i, create=False)

    def _link_network(self, create=False):
        s_set = set()
        e_set = set()
        r_set = set()
        for i in self.selfloop_df.to_dict("records"):
            s = ClusterNode(
                properties={
                    "cluster_index": i["clusterid1"],
                    "workflow_uuid": self.workflow_uuid,
                    "task": self.task,
                }
            )
            e = ClusterNode(
                properties={
                    "cluster_index": i["clusterid2"],
                    "workflow_uuid": self.workflow_uuid,
                    "task": self.task,
                }
            )
            r = MolecularNetwork(start=s, end=e)

            if "edgeannotation" not in i or i["edgeannotation"] == " ":
                i["edgeannotation"] = None
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

    def _link_library_to_chem(self):
        all_chem_nodes = set()
        all_rels = set()
        for i in self.library_nodes:
            chemstring = None
            if "inchi" in i.properties and i.properties["inchi"].startswith("InChI"):
                chemstring = i.properties["inchi"]
            elif "smiles" in i.properties and i.properties["smiles"]:
                chemstring = i.properties["smiles"]
            if chemstring:
                try:
                    cmpd = ChemicalCompound(chemstring)
                    a = cmpd.node
                    all_chem_nodes.add(a)
                    all_rels.add(GnpsLibraryToChem(start=i, end=a))
                except Exception as e:
                    log.debug(
                        f"Error creating chemical compound node for {i.properties['uid']}: {e}"
                    )
        all_chem_nodes = list(all_chem_nodes)
        all_rels = list(all_rels)
        if all_chem_nodes:
            all_chem_nodes[0].add_multiple_to_neo4j(all_chem_nodes, create=False)
        if all_rels:
            all_rels[0].add_multiple_to_neo4j(all_rels, create=False)
