import csv
import pickle
import tempfile
from multiprocessing import cpu_count
from operator import attrgetter
from pathlib import Path
from typing import List

from rich.progress import Progress

from socialgene.base.compare_protein import CompareProtein
from socialgene.base.molbio import Molbio
from socialgene.hmm.hmmer import HMMER
from socialgene.mmseqs.search import search as mmseqs_search
from socialgene.neo4j.neo4j import GraphDriver, Neo4jQuery
from socialgene.neo4j.search.basic import search_protein_hash
from socialgene.parsers.hmmer_parser import HmmerParser
from socialgene.parsers.sequence_parser import SequenceParser
from socialgene.utils.chunker import chunk_a_list_with_numpy
from socialgene.utils.file_handling import open_write
from socialgene.utils.logging import log


class SocialGene(Molbio, CompareProtein, SequenceParser, Neo4jQuery, HmmerParser):
    """Main class for building, sotring, working with protein and genomic info"""

    # _genomic_info_export_tablenames exports tables for the nextflow pipeline
    _genomic_info_export_tablenames = [
        "table_protein_to_go",
        "table_assembly_to_locus",
        "table_loci",
        "table_locus_to_protein",
        "table_assembly",
        "table_assembly_to_taxid",
        "table_protein_ids",
    ]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.protein_comparison = []

    ########################################
    # Get info
    ########################################

    def get_loci_ids(self):
        return list(self.loci.keys())

    ########################################
    # Filter
    ########################################

    def get_protein_domains_from_db(self, protein_id_list):
        # Search neo4j for a protein with the same hash, if it exists,
        # retrieve the domain info for it
        get_protein_domains_result = self.query_neo4j(
            cypher_name="get_protein_domains",
            param=protein_id_list,
        )
        temp = {
            "seq_pro_score": float(0),
            "evalue": float(0),
            "domain_bias": float(0),
            "domain_score": float(0),
            "seq_pro_bias": float(0),
            "hmm_from": int(0),
            "hmm_to": int(0),
            "ali_from": int(0),
            "ali_to": int(0),
            "env_from": int(0),
            "env_to": int(0),
        }
        if get_protein_domains_result:
            for protein in get_protein_domains_result:
                for domain in protein["domains"]:
                    new_dict = temp | domain["domain_properties"]
                    new_dict["hmm_id"] = domain["hmm_id"]
                    self.proteins[protein["p1.uid"]].add_domain(
                        **new_dict,
                    )
        else:
            log.info("None of the input protein ids were found in the database")

    def annotate_proteins_with_neo4j(
        self,
        protein_hash_ids: List[str] = None,
        chunk_size: int = 1000,
        annotate_all: bool = False,
        progress: bool = False,
    ):
        """
        The function `annotate_proteins_with_neo4j` takes a list of protein hash IDs, queries a database
        for matching proteins, and retrieves their HMM annotations.

        Args:
          protein_hash_ids (List[str]): A list of protein hash IDs. These are unique identifiers for proteins in the database that you want to annotate.
          chunk_size (int): The `chunk_size` parameter determines the number of proteins to query at a time. Proteins are divided into chunks to improve efficiency when querying the database.  Defaults to 1000
          annotate_all (bool): The `annotate_all` parameter is a boolean flag that determines whether to annotate all proteins in the database or not. If set to `True`, all proteins in the database will be annotated. If set to `False`, you need to provide a list of protein hash IDs in the `protein_hash_ids. Defaults to False
          progress (bool): The `progress` parameter is a boolean flag that determines whether or not to display a progress bar during the execution of the function. If `progress` is set to `True`, a progress bar will be displayed to track the progress of the function. If `progress` is set to `False`. Defaults to False

        """
        if protein_hash_ids is None and annotate_all:
            protein_hash_ids = self.get_all_protein_hashes()
        elif isinstance(protein_hash_ids, str):
            protein_hash_ids = [protein_hash_ids]
        search_result = search_protein_hash(protein_hash_ids)
        if not any([i for i in search_result.values()]):
            log.info("No identical proteins found in the database.")
        else:
            # log.info(
            #     f"{len([i for i in search_result.values() if i])} of {len(protein_hash_ids)} searched proteins were found in the database, pulling their HMM annotations into python..."
            # )
            # get the number of chunks needed to to query "chunk_size" proteins at a time
            prot_len = len(protein_hash_ids)
            n_chunks = prot_len // chunk_size
            if n_chunks == 0:
                n_chunks = 1
            if n_chunks > (len(protein_hash_ids) - 1):
                chunked_list = chunk_a_list_with_numpy(
                    input_list=protein_hash_ids, n_chunks=n_chunks
                )
            else:
                chunked_list = [protein_hash_ids]
            del protein_hash_ids
            if progress:
                with Progress(transient=True) as pg:
                    task = pg.add_task("Progress...", total=n_chunks)
                    for protein_id_list in chunked_list:
                        self.get_protein_domains_from_db(
                            protein_id_list=protein_id_list
                        )
                        pg.update(task, advance=1)
            else:
                for protein_id_list in chunked_list:
                    self.get_protein_domains_from_db(protein_id_list=protein_id_list)
        return search_result

    def search_with_mmseqs(self, target_database, argstring):
        with tempfile.NamedTemporaryFile() as temp_path:
            self.write_fasta(temp_path.name)
            self.mmseqs_results = mmseqs_search(
                fasta_path=temp_path.name,
                target_database=target_database,
                argstring=argstring,
            )

    def compare_to_another_sg_object(sg1, sg2, argstring=""):
        with tempfile.TemporaryDirectory() as tmpdirname:
            sg1fa_path = Path(tmpdirname, "sg1.faa")
            sg2fa_path = Path(tmpdirname, "sg2.faa")
            sg1.write_fasta(sg1fa_path)
            sg2.write_fasta(sg2fa_path)
            return mmseqs_search(
                fasta_path=str(sg1fa_path),
                target_database=str(sg2fa_path),
                argstring=argstring,
            )

    def annotate(
        self, use_neo4j_precalc: bool = False, neo4j_chunk_size: int = 1000, **kwargs
    ):
        """
        The `annotate` function is a convenience function that retrieves HMMER annotation for all
        proteins in a Socialgene object, either from a Neo4j database or by using HMMER directly.

        Args:
          use_neo4j_precalc (bool): A boolean flag indicating whether to use precalculated domain
        information from Neo4j for annotation. If set to True, the domains info will be retrieved from
        Neo4j, and proteins not found in Neo4j will be analyzed with HMMER. If set to False, all
        proteins will be analyzed with HMMER. Defaults to False
          neo4j_chunk_size (int): The `neo4j_chunk_size` parameter determines the number of proteins
        that will be sent to Neo4j at a time for annotation. Defaults to 1000
        """
        if use_neo4j_precalc:
            db_res = self.annotate_proteins_with_neo4j(
                annotate_all=True,
                chunk_size=neo4j_chunk_size,
            )
            n_not_in_db = [k for k, v in db_res.items() if not v]
        else:
            n_not_in_db = list(self.proteins.keys())
        if n_not_in_db:
            self.annotate_proteins_with_hmmscan(
                protein_id_list=n_not_in_db,
                **kwargs,
            )

    def annotate_proteins_with_hmmscan(
        self, protein_id_list, hmm_directory, cpus=None, **kwargs
    ):
        """Run hmmscan on Protein objects

        Args:
            protein_id_list (list): list of protein hash ids to run hmmscan on
            hmm_directory (str): path to directory containing hmm files
            cpus (int): number of parallel processes to spawn

        """
        log.info(f"Annotating {len(protein_id_list)} proteins with HMMER's hmmscan")
        if cpus is None:
            if cpu_count() > 2:
                cpus = cpu_count() - 1
            else:
                cpus = 1
        # create a list of proteins as FASTA
        temp1 = [
            self.single_protein_to_fasta(i)
            for i in protein_id_list
            if self.proteins[i].sequence
        ]
        if not temp1:
            log.info("None of the input sequences contain an amino acid sequence")
            return
        with tempfile.TemporaryDirectory() as tmpdirname:
            HMMER().hmmscan(
                fasta_path="-",
                input="\n".join(temp1).encode(),
                hmm_directory=hmm_directory,
                outdirectory=tmpdirname,
                overwrite=True,
                **kwargs,
            )
            # parse the resulting domtblout files, saving results to the class proteins/domains
            for i in Path(tmpdirname).glob("*.domtblout"):
                self.parse_hmmout(i, hmmsearch_or_hmmscan="hmmscan")

    def hydrate_from_proteins(self):
        """Given a SocialGene object with proteins, retrieve from a running Neo4j database all locus and asssembly info for those proteins"""
        for result in Neo4jQuery.query_neo4j(
            cypher_name="retrieve_protein_locations",
            param=list(self.proteins.keys()),
        ):
            self.add_assembly(result["assembly"])
            for locus in result["loci"]:
                _ = self.assemblies[result["assembly"]].add_locus(
                    external_id=locus["locus"]
                )
                for feature in locus["features"]:
                    _ = (
                        self.assemblies[result["assembly"]]
                        .loci[locus["locus"]]
                        .add_feature(
                            type="protein",
                            protein_hash=feature["protein_id"],
                            start=feature["locus_start"],
                            end=feature["locus_end"],
                            strand=feature["strand"],
                        )
                    )

    def fill_given_locus_range(self, locus_uid, start, end):
        """Given a locus uid that's in the database, pull asssembly, locus, protein info for those proteins"""

        with GraphDriver() as db:
            assembly_uid = db.run(
                """
                MATCH (a1:assembly)<-[:ASSEMBLES_TO]-(n1:nucleotide {uid: $locus_uid})
                RETURN a1.uid as a_uid
            """,
                locus_uid=locus_uid,
            ).single()
        if not assembly_uid:
            raise ValueError("No assembly found in database")
        else:
            assembly_uid = assembly_uid.value()
        self.add_assembly(assembly_uid)
        with GraphDriver() as db:
            res = db.run(
                """
                MATCH (n1:nucleotide {uid: $locus_uid})
                RETURN n1 as nucleotide_node
            """,
                locus_uid=locus_uid,
            ).single()
        if not res:
            raise ValueError(f"{locus_uid} not found in database")
        external_id = res.value().get("external_id")
        self.assemblies[assembly_uid].add_locus(external_id=external_id)

        self.assemblies[assembly_uid].loci[external_id].info = self.assemblies[
            assembly_uid
        ].loci[external_id].info | dict(res.value())

        with GraphDriver() as db:
            res = db.run(
                """
                MATCH (n1:nucleotide {uid: $locus_uid})-[e1:ENCODES]->(p1:protein)
                WHERE e1.start >= $start AND e1.start <= $end
                RETURN e1, p1
            """,
                locus_uid=locus_uid,
                start=start,
                end=end,
            )
            for feature in res:
                self.assemblies[assembly_uid].loci[external_id].add_feature(
                    type="protein",
                    protein_hash=feature.value().end_node["uid"],
                    **feature.value(),
                )
                _ = self.add_protein(hash_id=feature.value().end_node["uid"])
        # self.annotate_proteins_with_neo4j(
        #     protein_hash_ids=None, annotate_all=True, progress=False
        # )

    def hydrate_protein_info(self):
        """Pull name (original identifier) and description of proteins from Neo4j"""
        for result in Neo4jQuery.query_neo4j(
            cypher_name="get_protein_info",
            param=list(self.proteins.keys()),
        ):
            self.proteins[result["protein_id"]].external_protein_id = result["name"]
            self.proteins[result["protein_id"]].description = result["description"]

    ########################################
    # File Outputs
    ########################################

    def export_all_domains_as_tsv(self, outpath, **kwargs):
        """
        The function exports all domains as a TSV file, sorted by protein ID and mean envelope position.

        Args:
          outpath: The `outpath` parameter is the path to the output file where the TSV (Tab-Separated
        Values) data will be written. It specifies the location and name of the file.
        """
        _domain_counter = 0
        with open_write(outpath, **kwargs) as f:
            tsv_writer = csv.writer(f, delimiter="\t")
            # sort to standardize the write order
            ordered_prot_ids = list(self.proteins.keys())
            ordered_prot_ids.sort()
            for prot_id in ordered_prot_ids:
                # sort to standardize the write order
                for domain in self.proteins[
                    prot_id
                ].domain_list_sorted_by_mean_envelope_position:
                    _domain_counter += 1
                    _temp = [prot_id]
                    _temp.extend(list(domain.all_attributes.values()))
                    tsv_writer.writerow(_temp)
        log.info(f"Wrote {str(_domain_counter)} domains to {outpath}")

    def ferment_pickle(self, outpath):
        """
        The function `ferment_pickle` saves a SocialGene object to a Python pickle file.

        Args:
          outpath: The `outpath` parameter is a string that represents the path where the pickled object will be saved. It should include the file name and extension. For example, "/path/to/save/object.pickle".
        """
        with open(outpath, "wb") as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def eat_pickle(inpath):
        """
        The `eat_pickle` function reads a saved SocialGene pickle file from the given path and returns a SocialGene object.

        Args:
          inpath: The `inpath` parameter is a string that represents the path to the file from which we
        want to load the pickled object.

        Returns:
          SocialGene object
        """
        with open(inpath, "rb") as handle:
            return pickle.load(handle)

    def filter_proteins(self, hash_list: List):
        """Filter proteins by list of hash ids

        Args:
            hash_list (List): List of protein hash identifiers

        Returns:
            Generator: generator returning tuple of length two -> Generator[(str, 'Protein'])
        """
        return ((k, v) for k, v in self.proteins.items() if k in hash_list)

    def write_fasta(
        self,
        outpath,
        hash_list: List = None,
        external_protein_id: bool = False,
        **kwargs,
    ):
        """Write proteins to a FASTA file

        Args:
            outpath (str): path of file that FASTA entries will be appended to
            hash_list (List, optional): hash id of the protein(s) to export. Defaults to None.
            external_protein_id (bool, optional): Write protein identifiers as the hash (True) or the original identifier (False). Defaults to False.
            **kwargs: see print(open_write.__doc__)
        """

        with open_write(filepath=outpath, **kwargs) as handle:
            counter = 0
            if hash_list:
                temp_iter = self.filter_proteins(hash_list)
            else:
                temp_iter = self.proteins.items()
            for k, v in temp_iter:
                if v.sequence is not None:
                    counter += 1
                    if external_protein_id:
                        handle.writelines(v.fasta_string_defline_external_id)
                    else:
                        handle.writelines(v.fasta_string_defline_hash_id)

        log.info(f"Wrote {str(counter)} proteins to {outpath}")

    def single_protein_to_fasta(self, hash_id):
        """Create FASTA strings

        Args:
            hash_id (str): hash id of protein to export

        Returns:
            str: fasta string
        """
        return f">{self.proteins[hash_id].hash_id}\n{self.proteins[hash_id].sequence}"

    def write_n_fasta(self, outdir, n_splits=1, **kwargs):
        """
        The function `write_n_fasta` exports protein sequences split into multiple fasta files.

        Args:
          outdir: The `outdir` parameter is a string that specifies the directory where the fasta files
        will be saved.
          n_splits: The `n_splits` parameter in the `write_n_fasta` function determines the number of
        fasta files the protein sequences will be split into. By default, it is set to 1, meaning all
        the protein sequences will be written into a single fasta file. If you specify a value greater
        than. Defaults to 1
        """

        def split(a, n):
            # https://stackoverflow.com/a/2135920
            k, m = divmod(len(a), n)
            return (
                a[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)] for i in range(n)
            )

        protein_list = split(list(self.proteins.keys()), n_splits)
        counter = 1
        for protein_group in protein_list:
            with open_write(
                Path(outdir, f"fasta_split_{counter}.faa"), **kwargs
            ) as handle:
                for k, v in {
                    key: self.proteins.get(key) for key in sorted(protein_group)
                }.items():
                    handle.writelines(f">{k}\n{v.sequence}\n")
            counter += 1

    def _merge_proteins(self, sg_object):
        """
        The function merges the proteins from another SOcialGene object into the current object.

        Args:
          sg_object: The `sg_object` parameter is an object of the same class as the current object. It
        represents another object that contains a collection of proteins.
        """
        self.proteins.update(sg_object.proteins)

    def _merge_assemblies(self, sg_object):
        """
        The function `_merge_assemblies` merges assemblies and loci from `sg_object` into `self` if they
        do not already exist, and adds protein features to loci if they already exist.

        Args:
          sg_object: The `sg_object` parameter is an object of the same class as the current object.
        """
        for assembly_k, assembly_v in sg_object.assemblies.items():
            if assembly_k not in self.assemblies:
                self.assemblies[assembly_k] = assembly_v
            else:
                for locus_k, locus_v in assembly_v.loci.items():
                    if locus_k not in self.assemblies[assembly_k].loci:
                        self.assemblies[assembly_k].loci[locus_k] = locus_v
                    else:
                        for feature in (
                            sg_object.assemblies[assembly_k].loci[locus_k].features
                        ):
                            self.assemblies[assembly_k].loci[locus_k].add_feature(
                                type="protein",
                                protein_hash=feature.protein_hash,
                                start=feature.start,
                                end=feature.end,
                                strand=feature.strand,
                            )

    # def deep_merge_with_sg1(sg_object_1, sg_object_2):
    def __add__(self, sg_object):
        """
        The function merges proteins and assemblies from two objects and appends protein comparison data
        to a list or dataframe.

        Args:
          sg_object: The `sg_object` parameter is an object of the same class as the current object. It
        represents another instance of the class that you want to add to the current instance.
        """
        self._merge_proteins(sg_object)
        self._merge_assemblies(sg_object)
        self.protein_comparison.extend(sg_object.protein_comparison)

    # def break_loci(self, gap_tolerance: int = 10000):
    #     """
    #     The `break_loci` function breaks a locus (nucleotide sequence/contig/scaffold) into pieces based
    #     on a specified gap tolerance between proteins.

    #     Args:
    #       gap_tolerance (int): The `gap_tolerance` parameter is an optional integer that specifies the
    #     maximum allowed gap between proteins in nucleotides. If the gap between two proteins is greater
    #     than this value, the locus/contig/scaffold will be split into multiple parts. The default value
    #     for `gap_tolerance` is 100. Defaults to 10000
    #     """
    #     for assembly_key, assembly_value in self.assemblies.items():
    #         for loci_id in [i for i in assembly_value.loci.keys()]:
    #             assembly_value.loci[loci_id].sort_by_middle()
    #             # no need to split if there's only 1 feature
    #             if len(assembly_value.loci[loci_id].features) == 1:
    #                 log.info(
    #                     f"`len(loci_values.features) == 1`, doing nothing to {assembly_key}, {loci_id}"
    #                 )
    #             else:
    #                 current_int = 0
    #                 temp_locus = defaultdict(set)
    #                 # add the first feature
    #                 temp_locus[current_int].add(
    #                     assembly_value.loci[loci_id].features[0]
    #                 )
    #                 for j, k in zip(
    #                     assembly_value.loci[loci_id].features,
    #                     assembly_value.loci[loci_id].features[1:],
    #                 ):
    #                     if k.start - j.end > gap_tolerance:
    #                         current_int += 1
    #                     temp_locus[current_int].add(k)
    #                 if len(temp_locus) > 1:
    #                     for new_locus_k, new_locus_v in temp_locus.items():
    #                         new_name = f"{loci_id}_{new_locus_k}"
    #                         assembly_value.add_locus(new_name)
    #                         for i in new_locus_v:
    #                             assembly_value.loci[new_name].features.add(i)
    #                     del assembly_value.loci[loci_id]
    #                     log.info(
    #                         f"Broke assembly: {assembly_key} locus: {loci_id} into {len(temp_locus)} parts based on the provided gap_tolerance of {gap_tolerance}"
    #                     )
    #                     del temp_locus

    # # def prune_loci(self, min_genes=1):
    # #     """
    # #     The function `prune_loci` removes loci from an assembly if they have fewer features than the
    # #     specified minimum number of genes.

    # #     Args:
    # #       min_genes: The `min_genes` parameter is an optional parameter that specifies the minimum
    # #     number of genes a locus must have in order to be kept. If a locus has fewer genes than the
    # #     specified minimum, it will be removed from the `self.assemblies` dictionary. Defaults to 1
    # #     """
    # #     # remove loci if they don't have any features
    # #     assembly_keys = list(self.assemblies.keys())
    # #     for assembly_key in assembly_keys:
    # #         loci_keys = list(self.assemblies[assembly_key].loci.keys())
    # #         for loci_key in loci_keys:
    # #             if (
    # #                 len(self.assemblies[assembly_key].loci[loci_key].features)
    # #                 < min_genes
    # #             ):
    # #                 del self.assemblies[assembly_key].loci[loci_key]

    # def top_k_loci(self, topk=2):
    #     """
    #     The function `top_k_loci` removes loci from assemblies that are not in the top k based on the
    #     number of features they have.

    #     Args:
    #       topk: The parameter "topk" is an integer that specifies the number of top loci to keep in each
    #     assembly. Defaults to 2
    #     """
    #     assembly_keys = list(self.assemblies.keys())
    #     for assembly_key in assembly_keys:
    #         loci_dict = {
    #             k: len(v.features)
    #             for k, v in self.assemblies[assembly_key].loci.items()
    #         }
    #         to_delete = [
    #             i
    #             for i in list(self.assemblies[assembly_key].loci.keys())
    #             if i not in sorted(loci_dict, key=loci_dict.get, reverse=True)[:topk]
    #         ]
    #         for loci_key in to_delete:
    #             del self.assemblies[assembly_key].loci[loci_key]

    ########################################
    # OUTPUTS FOR NEXTFLOW
    ########################################

    def write_table(self, outdir: str, tablename: str, filename: str = None, **kwargs):
        """
        The function writes a table to a specified output directory in TSV format, for import into Neo4j.

        Args:
          outdir (str): The `outdir` parameter is a string that specifies the directory where the output
        file will be saved.
          tablename (str): The `tablename` parameter is a string that specifies the name of the table to
        be written.
          filename (str): The `filename` parameter is an optional argument that specifies the name of
        the file to be written. If no `filename` is provided, the `tablename` will be used as the
        filename.
        """
        if not filename:
            filename = tablename
        outpath = Path(outdir, filename)
        with open_write(outpath, **kwargs) as handle:
            tsv_writer = csv.writer(
                handle, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL
            )
            for i in getattr(self, tablename)():
                tsv_writer.writerow(i)

    def table_protein_to_go(self):
        """
        The function `table_protein_to_go` iterates through the assemblies, loci, and features of a
        given object, and yields tuples containing the protein hash and GO term for each feature that
        has GO terms, for import into Neo4j.
        """
        for av in self.assemblies.values():
            for v in av.loci.values():
                for feature in v.features:
                    # not all features will have goterms so check here
                    if feature.goterms:
                        for goterm in feature.goterms:
                            yield (
                                feature.protein_hash,
                                goterm.removeprefix("GO:").strip(),
                            )

    def table_locus_to_protein(self):
        """
        The function `table_locus_to_protein` generates a table of protein information from each locus, for import into Neo4j.
        """
        for ak, av in self.assemblies.items():
            for k, loci in av.loci.items():
                temp_list = list(loci.features)
                # sort features by id then start to maintain consistent output
                temp_list.sort(key=attrgetter("protein_hash"))
                temp_list.sort(key=attrgetter("start"))
                for feature in temp_list:
                    if feature.feature_is_protein():
                        yield (
                            feature.parent_object.uid,
                            feature.protein_hash,
                            feature.protein_id,
                            feature.locus_tag,
                            feature.start,
                            feature.end,
                            feature.strand,
                            feature.description,
                            feature.partial_on_complete_genome,
                            feature.missing_start,
                            feature.missing_stop,
                            feature.internal_stop,
                            feature.partial_in_the_middle_of_a_contig,
                            feature.missing_N_terminus,
                            feature.missing_C_terminus,
                            feature.frameshifted,
                            feature.too_short_partial_abutting_assembly_gap,
                            feature.incomplete,
                        )

    def table_protein_ids(self):
        """Protein hash id table for import into Neo4j"""
        for protein in self.proteins.values():
            yield (
                protein.hash_id,
                protein.crc64,
            )

    def table_assembly_to_locus(self):
        """Assembly to locus table for import into Neo4j

        Args:
            outdir (str, optional): Defaults to ".".
        """
        for ak, av in self.assemblies.items():
            for k, v in av.loci.items():
                #  ["assembly", "internal_locus_id"]
                yield (ak, v.uid)

    def table_loci(self):
        """
        Generate a table of loci information.

        Yields:
            tuple: (internal_locus_id, external_locus_id, [info])
        """
        for _av in self.assemblies.values():
            for locus in _av.loci.values():
                yield tuple(
                    [locus.uid]
                    + [locus.external_id]
                    + list(locus.metadata.all_attributes.values())
                )

    def table_assembly(self):
        """Assembly table for import into Neo4j

        Args:
            outdir (str, optional): Defaults to ".".
        """
        for assembly in self.assemblies.values():
            yield tuple(
                [assembly.uid] + list(assembly.metadata.all_attributes.values())
            )

    def table_assembly_to_taxid(self):
        """Assembly table for import into Neo4j

        Args:
            outdir (str, optional): Defaults to ".".
        """
        for k, v in self.assemblies.items():
            if v.taxid:
                yield (k, v.taxid)

    def drop_all_non_protein_features(self):
        """Drop features from all assembly/loci that aren't proteins/pseudo-proteins

        Returns:
            list: list of features that were removed (if return_removed=True)
        """
        for a_v in self.assemblies.values():
            for l_v in a_v.loci.values():
                l_v.drop_non_protein_features()
