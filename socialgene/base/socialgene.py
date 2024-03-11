import csv
import pickle
import tempfile
from multiprocessing import cpu_count
from operator import attrgetter
from pathlib import Path
from typing import List

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from rich.progress import Progress

from socialgene.base.compare_protein import CompareProtein
from socialgene.base.molbio import Molbio
from socialgene.hmm.hmmer import HMMER
from socialgene.neo4j.neo4j import GraphDriver, Neo4jQuery
from socialgene.neo4j.search.basic import search_protein_hash
from socialgene.parsers.hmmer_parser import HmmerParser
from socialgene.parsers.sequence_parser import SequenceParser
from socialgene.utils.chunker import chunk_a_list_with_numpy
from socialgene.utils.file_handling import open_write
from socialgene.utils.logging import log


class SocialGene(Molbio, CompareProtein, SequenceParser, Neo4jQuery, HmmerParser):
    """Main class for building, sotring, working with protein and genomic info"""

    # _export_table_names exports tables for the nextflow pipeline
    _export_table_names = [
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

    def get_all_gene_clusters(self):
        for i in self.assemblies.values():
            for j in i.loci.values():
                for k in j.gene_clusters:
                    yield k

    @property
    def protein_iter(self):
        for i in self.proteins.values():
            yield i

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
            log.info(
                "None of the input protein ids were found or had domains in the database"
            )

    def annotate_proteins_with_neo4j(
        self,
        protein_uids: List[str] = None,
        chunk_size: int = 1000,
        annotate_all: bool = False,
        progress: bool = False,
    ):
        """
        The function `annotate_proteins_with_neo4j` takes a list of protein hash IDs, queries a database
        for matching proteins, and retrieves their HMM annotations.

        Args:
          protein_uids (List[str]): A list of protein hash IDs. These are unique identifiers for proteins in the database that you want to annotate.
          chunk_size (int): The `chunk_size` parameter determines the number of proteins to query at a time. Proteins are divided into chunks to improve efficiency when querying the database.  Defaults to 1000
          annotate_all (bool): The `annotate_all` parameter is a boolean flag that determines whether to annotate all proteins in the database or not. If set to `True`, all proteins in the database will be annotated. If set to `False`, you need to provide a list of protein hash IDs in the `protein_uids. Defaults to False
          progress (bool): The `progress` parameter is a boolean flag that determines whether or not to display a progress bar during the execution of the function. If `progress` is set to `True`, a progress bar will be displayed to track the progress of the function. If `progress` is set to `False`. Defaults to False

        """
        if protein_uids is None and annotate_all:
            protein_uids = self.get_all_feature_uids()
        elif isinstance(protein_uids, str):
            protein_uids = [protein_uids]
        log.info(
            f"Searching database for HMM annotations of {len(protein_uids)} proteins."
        )
        search_result = search_protein_hash(protein_uids)
        if not any([i for i in search_result.values()]):
            log.info("No identical proteins found in the database.")
        else:
            # log.info(
            #     f"{len([i for i in search_result.values() if i])} of {len(protein_uids)} searched proteins were found in the database, pulling their HMM annotations into python..."
            # )
            # get the number of chunks needed to to query "chunk_size" proteins at a time
            prot_len = len(protein_uids)
            n_chunks = prot_len // chunk_size
            if n_chunks == 0:
                n_chunks = 1
            if n_chunks > (len(protein_uids) - 1):
                chunked_list = chunk_a_list_with_numpy(
                    input_list=protein_uids, n_chunks=n_chunks
                )
            else:
                chunked_list = [protein_uids]
            del protein_uids
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
        self, protein_id_list, hmm_filepath, cpus=None, **kwargs
    ):
        """Run hmmscan on Protein objects

        Args:
            protein_id_list (list): list of protein hash ids to run hmmscan on
            hmm_filepath (str): path to file of HMMs
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
            self.proteins[i].fasta_string_defline_uid
            for i in protein_id_list
            if self.proteins[i].sequence
        ]
        if not temp1:
            log.info("None of the input sequences contain an amino acid sequence")
            return
        with tempfile.NamedTemporaryFile() as temp_path:
            hmmer_temp = HMMER(hmm_filepath=hmm_filepath)
            hmmer_temp.hmmscan(
                fasta_path="-",
                input="\n".join(temp1).encode(),
                domtblout_path=temp_path.name,
                overwrite=True,
                cpus=cpus,
                **kwargs,
            )
            # parse the resulting domtblout file, saving results to the class proteins/domains
            self.parse_hmmout(temp_path.name, hmmsearch_or_hmmscan="hmmscan")

    def add_sequences_from_neo4j(self):
        try:
            with GraphDriver() as db:
                for i in db.run(
                    """
                    MATCH (p1:protein)
                    WHERE p1.uid in $uid
                    RETURN p1.uid as uid, p1.sequence as sequence
                    """,
                    uid=[i.uid for i in self.protein_iter],
                ):
                    if i["uid"] in self.proteins:
                        self.proteins[i["uid"]].sequence = i["sequence"]
        except Exception:
            log.debug(f"Error trying to retrieve domains for {self.uid}")

    def hydrate_from_proteins(self):
        """Given a SocialGene object with proteins, retrieve from a running Neo4j database all locus and asssembly info for those proteins"""
        for result in Neo4jQuery.query_neo4j(
            cypher_name="retrieve_protein_locations",
            param=list(self.proteins.keys()),
        ):
            self.add_assembly(uid=result["assembly"], parent=self)
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
                            uid=feature["external_id"],
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
        self.add_assembly(uid=assembly_uid, parent=self)
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
        self.assemblies[assembly_uid].loci[external_id].metadata.update(
            dict(res.value())
        )

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
                _ = self.add_protein(uid=feature.value().end_node["uid"])
                self.assemblies[assembly_uid].loci[external_id].add_feature(
                    type="protein",
                    uid=feature.value().end_node["uid"],
                    **feature.value(),
                )

        return {"assembly": assembly_uid, "locus": external_id}

    def _drop_all_cross_origin(self):
        # Drop likely cross-origin proteins (proteins that are longer than 100k bp)
        # have to remove from both gene cluster and locus
        for ak, av in self.assemblies.items():
            for nk, nv in av.loci.items():
                self.assemblies[ak].loci[nk].features = {
                    i for i in nv.features if abs(i.end - i.start) < 100000
                }
                for gene_cluster in self.assemblies[ak].loci[nk].gene_clusters:
                    gene_cluster.features = {
                        i
                        for i in gene_cluster.features
                        if abs(i.end - i.start) < 100000
                    }

    def hydrate_protein_info(self):
        """Pull name (original identifier) and description of proteins from Neo4j"""
        for result in Neo4jQuery.query_neo4j(
            cypher_name="get_protein_info",
            param=list(self.proteins.keys()),
        ):
            self.proteins[result["external_id"]].external_id = result["name"]
            self.proteins[result["external_id"]].description = result["description"]

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
                    _temp.extend(list(domain.all_attributes().values()))
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

    @property
    def fasta_string_defline_uid(self):
        for v in self.proteins.values():
            yield f">{v.uid}\n{v.sequence}\n"

    @property
    def fasta_string_defline_external_id(self):
        for v in self.proteins.values():
            yield f">{v.external_id}\n{v.sequence}\n"

    def write_fasta(
        self,
        outpath,
        external_id: bool = False,
        **kwargs,
    ):
        """Write all proteins to a FASTA file

        Args:
            outpath (str): path of file that FASTA entries will be appended to
            external_id (bool, optional): Write protein identifiers as the hash (True) or the original identifier (False). Defaults to False.
            **kwargs: see print(open_write.__doc__)
        """

        with open_write(filepath=outpath, **kwargs) as handle:
            counter = 0
            if external_id:
                fasta_gen = self.fasta_string_defline_external_id
            else:
                fasta_gen = self.fasta_string_defline_uid
            for i in fasta_gen:
                counter += 1
                handle.writelines(i)

        log.info(f"Wrote {str(counter)} proteins to {outpath}")

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

        # this can be done with itertools.batched in python 3.12
        def split(a, n):
            # https://stackoverflow.com/a/2135920
            k, m = divmod(len(a), n)
            return (
                a[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)] for i in range(n)
            )

        protein_list = split(
            [value for key, value in sorted(self.proteins.items(), reverse=False)],
            n_splits,
        )
        counter = 1
        for protein_group in protein_list:
            with open_write(
                Path(outdir, f"fasta_split_{counter}.faa"), **kwargs
            ) as handle:
                for i in protein_group:
                    handle.writelines(f">{i.uid}\n{i.sequence}\n")
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
                                uid=feature.uid,
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

    def write_genbank(self, outpath):
        for assembly in self.assemblies.values():
            for locus in assembly.loci.values():
                record = SeqRecord(
                    Seq(""),
                    id=locus.external_id,
                    name=locus.external_id,
                    description="A GenBank file generated by SocialGene.",
                    dbxrefs=[f"Assembly:{locus.parent.uid}"],
                )
                # Add annotation
                for feature in locus.features_sorted_by_midpoint:
                    biofeat = SeqFeature(
                        FeatureLocation(
                            start=feature.start,
                            end=feature.end,
                            strand=feature.strand,
                        ),
                        type=feature.type,
                        qualifiers={
                            k: v
                            for k, v in feature.all_attributes().items()
                            if v and k != "parent"
                        }
                        | {"translation": self.proteins[feature.uid].sequence},
                    )
                    record.features.append(biofeat)
                record.annotations["molecule_type"] = "DNA"

            SeqIO.write(
                record,
                outpath,
                "genbank",
            )

    def drop_all_non_protein_features(self, **kwargs):
        """Drop features from all assembly/loci that aren't proteins/pseudo-proteins

        Returns:
            list: list of features that were removed (if return_removed=True)
        """
        for a_v in self.assemblies.values():
            for l_v in a_v.loci.values():
                l_v.drop_non_protein_features()

    ########################################
    # OUTPUTS FOR NEXTFLOW
    ########################################

    def write_table(
        self,
        outdir: str,
        tablename: str,
        filename: str = None,
        include_sequences=False,
        **kwargs,
    ):
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
            for i in getattr(self, tablename)(include_sequences=include_sequences):
                tsv_writer.writerow(i)

    def table_protein_to_go(self, **kwargs):
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
                                feature.uid,
                                goterm.removeprefix("GO:").strip(),
                            )

    def table_locus_to_protein(self, **kwargs):
        """
        The function `table_locus_to_protein` generates a table of protein information from each locus, for import into Neo4j.
        """
        for ak, av in self.assemblies.items():
            for k, loci in av.loci.items():
                temp_list = list(loci.features)
                # sort features by id then start to maintain consistent output
                temp_list.sort(key=attrgetter("uid"))
                temp_list.sort(key=attrgetter("start"))
                for feature in temp_list:
                    if feature.feature_is_protein():
                        yield (
                            feature.parent.uid,
                            feature.uid,
                            feature.external_id,
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

    def table_protein_ids(self, include_sequences=False, **kwargs):
        """Protein hash id table for import into Neo4j"""
        for protein in self.proteins.values():
            if include_sequences:
                yield (protein.uid, protein.crc64, protein.sequence)
            else:
                yield (protein.uid, protein.crc64)

    def table_assembly_to_locus(self, **kwargs):
        """Assembly to locus table for import into Neo4j

        Args:
            outdir (str, optional): Defaults to ".".
        """
        for ak, av in self.assemblies.items():
            for k, v in av.loci.items():
                #  ["assembly", "internal_locus_id"]
                yield (ak, v.uid)

    def table_loci(self, **kwargs):
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
                    + list(locus.metadata.all_attributes().values())
                )

    def table_assembly(self, **kwargs):
        """Assembly table for import into Neo4j

        Args:
            outdir (str, optional): Defaults to ".".
        """
        for assembly in self.assemblies.values():
            yield tuple(
                [assembly.uid] + list(assembly.metadata.all_attributes().values())
            )

    def table_assembly_to_taxid(self, **kwargs):
        """Assembly table for import into Neo4j

        Args:
            outdir (str, optional): Defaults to ".".
        """
        for k, v in self.assemblies.items():
            if v.taxid:
                yield (k, v.taxid)
