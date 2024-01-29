from pathlib import Path
import json
import numpy as np
import pandas as pd
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log
from uuid import uuid4

import xml.etree.ElementTree as ET


def extract_parameter(xml_string):
    root = ET.fromstring(xml_string)
    return {root.attrib["name"]: root.text}


class GNPS_SNETS:
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
        self._filename_to_assembly()

    def _get_gnps_paths(self):
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
        # make a dict of the params ( lines look like: <parameter name="MAX_SHIFT">1999</parameter>)
        with open(self.params_xml_path, "r") as f:
            for i in f:
                try:
                    self.params.append(extract_parameter(i))
                except:
                    continue

    def _file_mapping(self):
        for i in self.params:
            for k, v in i.items():
                if k == "upload_file_mapping":
                    if (
                        v.split("|")[1].split("/")[0]
                        == [i["user"] for i in self.params if "user" in i][0]
                    ):
                        self.filemap[v.split("|")[0]] = v.split("|")[1].split("/")[-1]

    def _parse_clusterinfo(self):
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
        with open(self.params_xml_path, "r") as f:
            for i in f:
                if 'parameter name="uuid"' in i:
                    self.workflow_uuid = i.split(">")[1].split("<")[0]

    def _parse_dfs(self):
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

    def _filename_to_assembly(self):
        self.clusterinfo_df["assembly"] = self.clusterinfo_df[
            "original_filename"
        ].apply(lambda x: x.split("_")[0])

    def _check_db_for_assemblies(self):
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
        log.info("Adding GNPS library hits to db")
        with GraphDriver() as db:
            _ = db.run(
                """
                WITH $df as df
                UNWIND df as row
                MERGE (gls:gnps_library_spectrum {uid: row.spectrumid})
                    ON CREATE SET
                        gls.compound_name = row.compound_name,
                        gls.compound_source = row.compound_source,
                        gls.pi = row.pi,
                        gls.data_collector = row.data_collector,
                        gls.adduct = row.adduct,
                        gls.precursor_mz = row.precursor_mz,
                        gls.exactmass = row.exactmass,
                        gls.charge = row.charge,
                        gls.cas_number = row.cas_number,
                        gls.pubmed_id = row.pubmed_id,
                        gls.smiles = row.smiles,
                        gls.inchi = row.inchi,
                        gls.inchi_aux = row.inchi_aux,
                        gls.library_class = row.library_class,
                        gls.ionmode = row.ionmode,
                        gls.libraryqualitystring = row.libraryqualitystring,
                        gls.mqscore = row.mqscore,
                        gls.tic_query = row.tic_query,
                        gls.rt_query = row.rt_query,
                        gls.mzerrorppm = row.mzerrorppm,
                        gls.sharedpeaks = row.sharedpeaks,
                        gls.massdiff = row.massdiff,
                        gls.libmz = row.libmz,
                        gls.specmz = row.specmz,
                        gls.speccharge = row.speccharge,
                        gls.moleculeexplorerdatasets = row.moleculeexplorerdatasets,
                        gls.moleculeexplorerfiles = row.moleculeexplorerfiles,
                        gls.molecular_formula = row.molecular_formula,
                        gls.inchikey = row.inchikey,
                        gls.inchikey_planar = row.inchikey_planar
                """,
                df=self.specnets_df.to_dict("records"),
            ).value()

    def add_gnps_library_spectrum_classifications(self):
        log.info("Adding GNPS library classifications to db")
        with GraphDriver() as db:
            results = db.run(
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
        log.info("Adding GNPS library ion sources to db")
        with GraphDriver() as db:
            results = db.run(
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
        log.info("Adding GNPS library instruments to db")
        with GraphDriver() as db:
            results = db.run(
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
        log.info("Adding GNPS library organisms to db")
        with GraphDriver() as db:
            results = db.run(
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
        log.info("Adding GNPS clusters to db")
        with GraphDriver() as db:
            results = db.run(
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
        log.info("Adding links between GNPS clusters and assemblies")
        df = (
            self.clusterinfo_df.groupby(["assembly"])["clusteridx"]
            .apply(set)
            .apply(list)
        )
        for i in df.reset_index().itertuples(index=False):
            with GraphDriver() as db:
                results = db.run(
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
        log.info("Adding and linking input spectra to GNPS clusters")
        with GraphDriver() as db:
            results = db.run(
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
        log.info("Adding links between GNPS clusters")
        with GraphDriver() as db:
            results = db.run(
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
