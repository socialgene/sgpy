from socialgene.utils.logging import log
from pathlib import Path
import csv

# `SocialgeneModules` groups `Neo4jImportData`
# `Neo4jImportData` defines information and structure about data that will be imported into Neo4j


hmm_sources = [
    "pfam",
    "antismash",
    "tigrfam",
    "amrfinder",
    "prism",
    "resfams",
    "bigslice",
    "classiphage",
    "virus_orthologous_groups",
    "local",
]


def parse_hmmlist_input(input):
    if input == "all" or "all" in input:
        temp = [i for i in hmm_sources if i != "local"]
        return temp
    else:
        temp = [i for i in input if i in hmm_sources]
        return temp


class Neo4jImportData:
    def __init__(
        self,
    ):
        self.nodes = {
            "parameters": {
                "neo4j_label": "parameters",
                "header_filename": "parameters.header",
                "target_subdirectory": "parameters",
                "target_extension": "socialgene_parameters",
                "header": [
                    "id:ID(when)",
                    "SG_LOC_NEO4J",
                    "SG_LOC_HMMS",
                    "NEO4J_dbms_memory_pagecache_size",
                    "NEO4J_dbms_memory_heap_initial__size",
                    "NEO4J_dbms_memory_heap_max__size",
                    "HMMSEARCH_IEVALUE",
                    "HMMSEARCH_BACKGROUND",
                    "HMMSEARCH_BIASFILTER",
                    "HMMSEARCH_NULL2",
                    "HMMSEARCH_SEED",
                    "HMMSEARCH_Z",
                    "HMMSEARCH_DOMZ",
                    "HMMSEARCH_F1",
                    "HMMSEARCH_F2",
                    "HMMSEARCH_F3",
                    "HMMSEARCH_E",
                    "HMMSEARCH_DOME",
                    "HMMSEARCH_INCE",
                    "HMMSEARCH_INCDOME",
                    "HMMSEARCH_BITCUTOFFS",
                    "platform",
                    "architecture",
                    "py_executable",
                    "py_version",
                    "genome_download_command",
                ],
            },
            "assembly": {
                "neo4j_label": "assembly",
                "header_filename": "assembly.header",
                "target_subdirectory": "genomic_info",
                "target_extension": "assemblies",
                "header": [
                    "id:ID(assembly)",
                    "mol_type",
                    "altitude",
                    "bio_material",
                    "cell_line",
                    "cell_type",
                    "chromosome",
                    "clone",
                    "clone_lib",
                    "collected_by",
                    "collection_date",
                    "country",
                    "cultivar",
                    "culture_collection",
                    "db_xref",
                    "dev_stage",
                    "ecotype",
                    "environmental_sample",
                    "focus",
                    "germline",
                    "haplogroup",
                    "haplotype",
                    "host",
                    "identified_by",
                    "isolate",
                    "isolation_source",
                    "lab_host",
                    "lat_lon",
                    "macronuclear",
                    "map",
                    "mating_type",
                    "metagenome_source",
                    "note",
                    "organelle",
                    "PCR_primers",
                    "plasmid",
                    "pop_variant",
                    "proviral",
                    "rearranged",
                    "segment",
                    "serotype",
                    "serovar",
                    "sex",
                    "specimen_voucher",
                    "strain",
                    "sub_clone",
                    "submitter_seqid",
                    "sub_species",
                    "sub_strain",
                    "tissue_lib",
                    "tissue_type",
                    "transgenic",
                    "type_material",
                    "variety",
                ],
            },
            "nucleotide": {
                "neo4j_label": "nucleotide",
                "header_filename": "locus.header",
                "target_subdirectory": "genomic_info",
                "target_extension": "loci",
                "header": [
                    "internal_id:ID(nucleotide)",
                    "id",
                    "mol_type",
                    "altitude",
                    "bio_material",
                    "cell_line",
                    "cell_type",
                    "chromosome",
                    "clone",
                    "clone_lib",
                    "collected_by",
                    "collection_date",
                    "country",
                    "cultivar",
                    "culture_collection",
                    "db_xref",
                    "dev_stage",
                    "ecotype",
                    "environmental_sample",
                    "focus",
                    "germline",
                    "haplogroup",
                    "haplotype",
                    "host",
                    "identified_by",
                    "isolate",
                    "isolation_source",
                    "lab_host",
                    "lat_lon",
                    "macronuclear",
                    "map",
                    "mating_type",
                    "metagenome_source",
                    "note",
                    "organelle",
                    "PCR_primers",
                    "plasmid",
                    "pop_variant",
                    "proviral",
                    "rearranged",
                    "segment",
                    "serotype",
                    "serovar",
                    "sex",
                    "specimen_voucher",
                    "strain",
                    "sub_clone",
                    "submitter_seqid",
                    "sub_species",
                    "sub_strain",
                    "tissue_lib",
                    "tissue_type",
                    "transgenic",
                    "type_material",
                    "variety",
                ],
            },
            "protein": {
                "neo4j_label": "protein",
                "header_filename": "protein_info.header",
                "target_subdirectory": "protein_info",
                "target_extension": "protein_info",
                "header": [
                    "id:ID(protein)",
                    # "source_db",
                    "name",
                    "description",
                    "seqlen:int",
                ],
            },
            "hmm": {
                "neo4j_label": "hmm",
                "header_filename": "sg_hmm_nodes_out.header",
                "target_subdirectory": "hmm_tsv_parse",
                "target_extension": "sg_hmm_nodes_out",
                "header": ["id:ID(hmm)", "model_length"],
            },
            "goterm": {
                "neo4j_label": "goterm",
                "header_filename": "goterm.header",
                "target_subdirectory": "tigrfam_info",
                "target_extension": "goterm",
                "header": [
                    "id:ID(goterm)",
                ],
            },
            "tigrfam_role": {
                "neo4j_label": "tigrfam_role",
                "header_filename": "tigrfam_role.header",
                "target_subdirectory": "tigrfam_info",
                "target_extension": "tigrfam_role",
                "header": ["id:ID(tigrfam_role)"],
            },
            "tigrfam_mainrole": {
                "neo4j_label": "tigrfam_mainrole",
                "header_filename": "tigrfam_mainrole.header",
                "target_subdirectory": "tigrfam_info",
                "target_extension": "tigrfam_mainrole",
                "header": ["id:ID(tigrfam_mainrole)"],
            },
            "tigrfam_subrole": {
                "neo4j_label": "tigrfam_subrole",
                "header_filename": "tigrfam_subrole.header",
                "target_subdirectory": "tigrfam_info",
                "target_extension": "tigrfam_subrole",
                "header": ["id:ID(tigrfam_subrole)"],
            },
            "taxid": {
                "neo4j_label": "taxid",
                "header_filename": "taxid.header",
                "target_subdirectory": "taxdump_process",
                "target_extension": "nodes_taxid",
                "header": ["id:ID(taxid)", "name", "rank"],
            },
            "mz_cluster_index": {
                "neo4j_label": "mz_cluster_index",
                "header_filename": "mz_cluster_index_nodes.header",
                "target_subdirectory": "paired_omics",
                "target_extension": "mz_cluster_index_nodes",
                "header": [
                    "id:ID(mz_cluster_index)",
                    "component_index:int",
                    "parent_mass:double",
                    "precursor_mass:double",
                    "sum_precursor_intensity:double",
                    "Smiles:String",
                    "rt_mean::double",
                    "rt_std_err:double",
                    "library_id:String",
                    "mq_score:double",
                    "mz_error_ppm:double",
                    "mass_diff:String",
                ],
            },
            "mz_source_file": {
                "neo4j_label": "mz_source_file",
                "header_filename": "mz_source_file.header",
                "target_subdirectory": "paired_omics",
                "target_extension": "hash.mz_source_file",
                "header": [
                    "id:ID(mz_source_file)",
                ],
            },
        }
        self.relationships = {
            "annotates": {
                "neo4j_label": "ANNOTATES",
                "header_filename": "protein_to_hmm_header.header",
                "target_subdirectory": "parsed_domtblout",
                "target_extension": "parseddomtblout",
                "header": [
                    ":END_ID(protein)",
                    ":START_ID(hmm)",
                    "env_from:int",
                    "env_to:int",
                    "seq_pro_score:float",
                    "evalue:float",
                    "i_evalue:float",
                    "domain_bias:float",
                    "domain_score:float",
                    "seq_pro_bias:float",
                    "hmm_from:int",
                    "hmm_to:int",
                    "ali_from:int",
                    "ali_to:int",
                ],
            },
            "assembles_to": {
                "neo4j_label": "ASSEMBLES_TO",
                "header_filename": "assembly_to_locus.header",
                "target_subdirectory": "genomic_info",
                "target_extension": "assembly_to_locus",
                "header": [":END_ID(assembly)", ":START_ID(nucleotide)"],
            },
            "contains": {
                "neo4j_label": "ENCODES",
                "header_filename": "locus_to_protein.header",
                "target_subdirectory": "genomic_info",
                "target_extension": "locus_to_protein",
                "header": [
                    ":START_ID(nucleotide)",
                    ":END_ID(protein)",
                    "start:int",
                    "end:int",
                    "strand:int",
                ],
            },
            "belongs_to": {
                "neo4j_label": "BELONGS_TO",
                "header_filename": "taxid_to_taxid.header",
                "target_subdirectory": "taxdump_process",
                "target_extension": "taxid_to_taxid",
                "header": [":START_ID(taxid)", ":END_ID(taxid)"],
            },
            "go_ann": {
                "neo4j_label": "GO_ANN",
                "header_filename": "tigrfam_to_go.header",
                "target_subdirectory": "tigrfam_info",
                "target_extension": "tigrfam_to_go",
                "header": [":START_ID(tigrfam)", ":END_ID(goterm)"],
            },
            "role_ann": {
                "neo4j_label": "ROLE_ANN",
                "header_filename": "tigrfam_to_role.header",
                "target_subdirectory": "tigrfam_info",
                "target_extension": "tigrfam_to_role",
                "header": [":START_ID(tigrfam)", ":END_ID(tigrfam_role)"],
            },
            "mainrole_ann": {
                "neo4j_label": "MAINROLE_ANN",
                "header_filename": "tigrfamrole_to_mainrole.header",
                "target_subdirectory": "tigrfam_info",
                "target_extension": "tigrfamrole_to_mainrole",
                "header": [":START_ID(tigrfam_role)", ":END_ID(tigrfam_mainrole)"],
            },
            "subrole_ann": {
                "neo4j_label": "SUBROLE_ANN",
                "header_filename": "tigrfamrole_to_subrole.header",
                "target_subdirectory": "tigrfam_info",
                "target_extension": "tigrfamrole_to_subrole",
                "header": [":START_ID(tigrfam_role)", ":END_ID(tigrfam_subrole)"],
            },
            "assembly_to_taxid": {
                "neo4j_label": "TAXONOMY",
                "header_filename": "assembly_to_taxid.header",
                "target_subdirectory": "genomic_info",
                "target_extension": "assembly_to_taxid",
                "header": [":START_ID(assembly)", ":END_ID(taxid)"],
            },
            "blastp": {
                "neo4j_label": "BLASTP",
                "header_filename": "blastp.header",
                "target_subdirectory": "diamond_blastp",
                "target_extension": "blast6",
                "header": [
                    ":START_ID(protein)",
                    ":END_ID(protein)",
                    "pident:float",
                    "length:int",
                    "mismatch:int",
                    "gapopen:int",
                    "qstart:int",
                    "qend:int",
                    "sstart:int",
                    "send:int",
                    "evalue:float",
                    "bitscore:float",
                    "qcovhsp:float",
                ],
            },
            "mmseqs2": {
                "neo4j_label": "MMSEQS2",
                "header_filename": "mmseqs2.header",
                "target_subdirectory": "mmseqs2_easycluster",
                "target_extension": "mmseqs2_results_cluster.tsv",
                "header": [":START_ID(protein)", ":END_ID(protein)"],
            },
            "cluster_to_file": {
                "neo4j_label": "CLUSTER_TO_FILE",
                "header_filename": "cluster_to_source_file.header",
                "target_subdirectory": "paired_omics",
                "target_extension": "cluster_to_source_file",
                "header": [":END_ID(mz_cluster_index)", ":START_ID(mz_source_file)"],
            },
            "molecular_network": {
                "neo4j_label": "MOLECULAR_NETWORK",
                "header_filename": "molecular_network.header",
                "target_subdirectory": "paired_omics",
                "target_extension": "molecular_network",
                "header": [
                    ":START_ID(mz_cluster_index)",
                    ":END_ID(mz_cluster_index)",
                    "delta_mz:double",
                    "meh:float",
                    "cosine:float",
                    "other_score:float",
                ],
            },
            "metabo": {
                "neo4j_label": "METABO",
                "header_filename": "assembly_to_mz_file.header",
                "target_subdirectory": "paired_omics",
                "target_extension": "assembly_to_mz_file",
                "header": [":START_ID(assembly)", ":END_ID(mz_source_file)"],
            },
        }
        # create hmm source node/relationship headers
        for i in hmm_sources:
            self.nodes[i] = {
                "neo4j_label": i,
                "header_filename": f"{i}_hmms_out.header",
                "target_subdirectory": "hmm_tsv_parse",
                "target_extension": f"{i}_hmms_out",
                "header": [
                    ":IGNORE",
                    "accession",
                    f"id:ID({i})",
                    "description",
                    "category",
                ],
            }
            self.relationships[i] = {
                "neo4j_label": "SOURCE_DB",
                "header_filename": f"{i}_hmms_out_relationships.header",
                "target_subdirectory": "hmm_tsv_parse",
                "target_extension": f"{i}_hmms_out",
                "header": [
                    ":START_ID(hmm)",
                    ":IGNORE",
                    f":END_ID({i})",
                    ":IGNORE",
                    ":IGNORE",
                ],
            }


class SocialgeneModules:
    def __init__(
        self,
    ):
        # if adding both a node and relationship for the same module, give it the same name in both dicts
        self.nodes = {
            "base": [
                "parameters",
                "assembly",
                "nucleotide",
                "protein",
            ],
            "hmms": [],
            "ncbi_taxonomy": ["taxid"],
            "base_hmm": ["hmm"],
            "tigrfam": [
                "goterm",
                "tigrfam_mainrole",
                "tigrfam_subrole",
                "tigrfam_role",
            ],
            "paired_omics": ["mz_cluster_index", "mz_source_file"],
        }
        self.relationships = {
            "base": [
                "contains",
                "assembles_to",
            ],
            "base_hmm": ["annotates"],
            "hmms": [],
            "tigrfam": [
                "mainrole_ann",
                "role_ann",
                "subrole_ann",
                "go_ann",
            ],
            "ncbi_taxonomy": ["belongs_to", "assembly_to_taxid"],
            "paired_omics": ["cluster_to_file", "molecular_network", "metabo"],
        }
        self.nodes["hmms"].extend(hmm_sources)
        self.relationships["hmms"].extend(hmm_sources)
        # enable all node and relationships as individual sgmodules options
        for i in Neo4jImportData().nodes.keys():
            if i not in self.nodes:
                self.nodes[i] = [i]
        for i in Neo4jImportData().relationships.keys():
            if i not in self.relationships:
                self.relationships[i] = [i]
        self.sg_modules = {"nodes": self.nodes, "relationships": self.relationships}

    def node_keylist(self):
        return list(self.nodes.keys())

    def relationship_keylist(self):
        return list(self.relationships.keys())

    @staticmethod
    def _filter(input_sg_modules, input_dict):
        temp = {k: v for k, v in input_dict.items() if k in input_sg_modules}
        # had warning for user if node/relationship is missing, but can't since not all groups have a node or vice versa
        return temp

    def filter_nodes(self, input_sg_modules):
        if "hmms" in input_sg_modules:
            input_sg_modules.append("base_hmm")
        return self._filter(input_sg_modules, self.nodes)

    def filter_relationships(self, input_sg_modules):
        if "hmms" in input_sg_modules:
            input_sg_modules.append("base_hmm")
        return self._filter(input_sg_modules, self.relationships)


def _writer(outdir, header_dict):
    outpath = Path(outdir, f"{header_dict['header_filename']}")
    with open(outpath, "w") as tsv_output_con:
        tsv_writer = csv.writer(
            tsv_output_con,
            delimiter="\t",
            quotechar='"',
            quoting=csv.QUOTE_MINIMAL,
        )
        tsv_writer.writerow(header_dict["header"])
        log.info(f"\tWriting {header_dict['header_filename']} to: {outpath}")


def write_neo4j_headers(sg_modules: list, hmmlist: list, outdir: str):
    sg_mod_object = SocialgeneModules()
    header_object = Neo4jImportData()
    reduced_dict = {
        "nodes": sg_mod_object.filter_nodes(sg_modules),
        "relationships": sg_mod_object.filter_relationships(sg_modules),
    }
    for rd_key, rd_value in reduced_dict.items():
        for sg_mod_key, header_keys in rd_value.items():
            if sg_mod_key == "hmms":
                keylist = [i for i in header_keys if i in hmmlist]
                keylist.extend(rd_value["base_hmm"])
            else:
                keylist = header_keys
            for i in keylist:
                _writer(outdir, getattr(header_object, rd_key)[i])
