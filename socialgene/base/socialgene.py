# python dependencies
import pickle
from pathlib import Path
import csv
from copy import deepcopy
import tempfile
from multiprocessing import cpu_count
import itertools
from functools import partial
import gzip
from collections import defaultdict
from operator import attrgetter

# external dependencies
from rich.progress import Progress
import pandas as pd
from socialgene.hmm.hmmer import HMMER

# internal dependencies
from socialgene.parsers.sequence_parser import SequenceParser
from socialgene.parsers.hmmer_parser import HmmerParser
from socialgene.base.molbio import Molbio
from socialgene.utils.logging import log
from socialgene.base.compare_protein import CompareProtein

from socialgene.neo4j.neo4j import Neo4jQuery

from socialgene.scoring.scoring import mod_score
import socialgene.hashing.hashing as hasher
from socialgene.utils.chunker import chunk_a_list_with_numpy
from socialgene.neo4j.search.basic import search_protein_hash


class SocialGene(Molbio, CompareProtein, SequenceParser, Neo4jQuery, HmmerParser):
    """Main class for building, sotring, working with protein and genomic info"""

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
        self.protein_comparison = []

    ########################################
    # Get info
    ########################################

    def get_loci_ids(self):
        return list(self.loci.keys())

    def get_min_maxcoordinates(self):
        return {k: v.get_min_maxcoordinates() for k, v in self.assemblies.items()}

    ########################################
    # Filter
    ########################################
    def filter_domains(self):
        """Prune all domains in all proteins that don't meet the inclusion threshold (currently HMMER's i_evalue)"""
        _before_count = sum([len(i.domains) for i in self.proteins.values()])
        for protein in self.proteins.values():
            for domain in protein.domains:
                if not (domain.domain_within_threshold()):
                    protein.domains.remove(domain)
        _after_count = sum([len(i.domains) for i in self.proteins.values()])
        log.info(
            f"Removed {str(_before_count - _after_count)} domains. Before: {str(_before_count)}; After: {str(_after_count)}"
        )

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
                    self.proteins[protein["p1.id"]].add_domain(
                        exponentialized=False,
                        **new_dict,
                    )
        else:
            log.info("None of the input protein ids were found in the database")

    def annotate_proteins_with_neo4j(
        self,
        protein_hash_ids=None,
        chunk_size=1000,
        annotate_all=False,
        progress=False,
    ):

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

    def annotate(
        self, use_neo4j_precalc: bool = False, neo4j_chunk_size: int = 1000, **kwargs
    ):
        """Convenience function to get HMMER annotation for all proteins in a Socialgene() object

        Args:
            use_neo4j_precalc (bool, optional): If True, domains info will be retrieved from Neo4j, proteins not in Neo4j will be analyzed with HMMER. Defaults to False.
            neo4j_chunk_size (int, optional): Send x-proteins at a time to Neo4j for annotation
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
            self.annotate_with_hmmscan(
                protein_id_list=n_not_in_db,
                **kwargs,
            )

    def annotate_with_hmmscan(self, protein_id_list, hmm_filepath, cpus=None, **kwargs):
        """Run hmmscan in parallel (meant to be when the number of proteins is relatively small)

        Args:
            protein_id_list (list): list of protein hash ids to run hmmscan on
            hmm_filepath (str): path to hmm file (must have be hmmpressed first)
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
        # create a list of lists of proteins as FASTA
        if cpus > (len(temp1) - 1) and len(temp1) > 10:
            temp2 = chunk_a_list_with_numpy(input_list=temp1, n_chunks=cpus)
        else:
            temp2 = [temp1]
        del temp1
        # create tempfiles for hmmscan to write domtblout files to
        files = [tempfile.NamedTemporaryFile() for i in temp2]
        filenames = [i.name for i in files]
        for i in zip(temp2, itertools.repeat(hmm_filepath), filenames):
            HMMER.hmmscan(
                fasta_path="-",
                input="\n".join(i[0]).encode(),
                hmm_filepath=i[1],
                domtblout_path=i[2],
                overwrite=True,
                **kwargs,
            )
        # parse the resulting domtblout files, saving results to the class proteins/domains
        for i, ii in zip(filenames, files):
            print(i)
            self.parse_hmmout(i, hmmsearch_or_hmmscan="hmmscan")
            ii.close()

    @staticmethod
    def compare_two_proteins(protein_1, protein_2):
        return mod_score(protein_1.get_domain_vector(), protein_2.get_domain_vector())

    def fill_locus_assembly_from_db(self):
        """Given a SocialGene object with proteins, retrieve from a running Neo4j database all locus and asssembly info for those proteins"""
        for result in Neo4jQuery.query_neo4j(
            cypher_name="retrieve_protein_locations",
            param=list(self.proteins.keys()),
        ):
            self.add_assembly(result["assembly"])
            for locus in result["loci"]:
                _ = self.assemblies[result["assembly"]].add_locus(id=locus["locus"])
                for feature in locus["features"]:
                    _ = (
                        self.assemblies[result["assembly"]]
                        .loci[locus["locus"]]
                        .add_feature(
                            type="protein",
                            id=feature["protein_id"],
                            start=feature["locus_start"],
                            end=feature["locus_end"],
                            strand=feature["strand"],
                        )
                    )

    def get_db_protein_info(self):
        """Pull name (original identifier) and description of proteins from Neo4j"""
        for result in Neo4jQuery.query_neo4j(
            cypher_name="get_protein_info",
            param=list(self.proteins.keys()),
        ):
            self.proteins[result["protein_id"]].other_id = result["name"]
            self.proteins[result["protein_id"]].description = result["description"]

    ########################################
    # File Outputs
    ########################################

    def export_all_domains_as_tsv(self, outpath):
        _domain_counter = 0
        with open(outpath, "a") as f:
            tsv_writer = csv.writer(f, delimiter="\t")
            # sort to standardize the write order
            ordered_prot_ids = list(self.proteins.keys())
            ordered_prot_ids.sort()
            for prot_id in ordered_prot_ids:
                # sort to standardize the write order
                for domain in self.proteins[
                    prot_id
                ].sort_domains_by_mean_envelope_position():
                    _domain_counter += 1
                    _temp = [prot_id]
                    _temp.extend(list(domain.get_dict().values()))
                    tsv_writer.writerow(_temp)
        log.info(f"Wrote {str(_domain_counter)} domains to {outpath}")

    def ferment_pickle(self, outpath):
        with open(outpath, "wb") as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def eat_pickle(inpath):
        with open(inpath, "rb") as handle:
            return pickle.load(handle)

    def write_fasta(self, outpath, other_id=False, gz=False):
        if gz:
            _open = partial(gzip.open, mode="at")
        else:
            _open = partial(open, mode="a")
        with _open(outpath) as handle:
            counter = 0
            for k, v in self.proteins.items():
                if v.sequence is not None:
                    counter += 1
                    if other_id:
                        out_id = self.proteins[k].other_id
                    else:
                        out_id = k
                    handle.writelines(f">{out_id}\n{v.sequence}\n")
        log.info(f"Wrote {str(counter)} proteins to {outpath}")

    def single_protein_to_fasta(self, hash_id):
        """Create FASTA strings

        Args:
            hash_id (str): hash id of protein to export

        Returns:
            str: fasta string
        """
        return f">{self.proteins[hash_id].hash_id}\n{self.proteins[hash_id].sequence}"

    def write_n_fasta(self, outdir, n_splits=1):
        """Export protein sequences split into n-number of fasta files

        Args:
            outdir (str): Directory to save fasta files into
            n_splits (int, optional): Defaults to 1.
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
            with open(Path(outdir, f"fasta_split_{counter}.faa"), "w") as handle:
                for k, v in {
                    key: self.proteins.get(key) for key in protein_group
                }.items():
                    handle.writelines(f">{k}\n{v.sequence}\n")
            counter += 1

    def _merge_proteins(self, sg_object):
        self.proteins.update(sg_object.proteins)

    def _merge_assemblies(self, sg_object):
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
                                id=feature.id,
                                start=feature.start,
                                end=feature.end,
                                strand=feature.strand,
                            )

    # def deep_merge_with_sg1(sg_object_1, sg_object_2):
    def __add__(self, sg_object):
        self._merge_proteins(sg_object)
        self._merge_assemblies(sg_object)
        # check if self.protein_comparison is pandas or a list and append accordingly
        if isinstance(self.protein_comparison, pd.core.frame.DataFrame):
            sg_object.protein_comparison_to_df()
            self.protein_comparison.append(sg_object.protein_comparison)
        elif isinstance(self.protein_comparison, list):
            self.protein_comparison.extend(sg_object.protein_comparison)
        else:
            raise TypeError()

    def break_loci(self, gap_tolerance: int = 10000):
        """Break a locus (nucleotide sequence/contig/scaffold) into pieces based on gap tolerance betwen proteins (measured in nucleotides not amino acids)
        Args:
            gap_tolerance (int, optional): Allowed nucleotide gap between proteins, any greater and the locus/contig/scaffold will be split. Defaults to 10000.
        """
        for assembly_key, assembly_value in self.assemblies.items():
            for loci_id in [i for i in assembly_value.loci.keys()]:
                assembly_value.loci[loci_id].sort_by_middle()
                # no need to split if there's only 1 feature
                if len(assembly_value.loci[loci_id].features) == 1:
                    log.info(
                        f"`len(loci_values.features) == 1`, doing nothing to {assembly_key}, {loci_id}"
                    )
                else:
                    current_int = 0
                    temp_locus = defaultdict(set)
                    # add the first feature
                    temp_locus[current_int].add(
                        assembly_value.loci[loci_id].features[0]
                    )
                    for j, k in zip(
                        assembly_value.loci[loci_id].features,
                        assembly_value.loci[loci_id].features[1:],
                    ):
                        if k.start - j.end > gap_tolerance:
                            current_int += 1
                        temp_locus[current_int].add(k)
                    if len(temp_locus) > 1:
                        for new_locus_k, new_locus_v in temp_locus.items():
                            new_name = f"{loci_id}_{new_locus_k}"
                            assembly_value.add_locus(new_name)
                            for i in new_locus_v:
                                assembly_value.loci[new_name].features.add(i)
                        del assembly_value.loci[loci_id]
                        log.info(
                            f"Broke assembly: {assembly_key} locus: {loci_id} into {len(temp_locus)} parts based on the provided gap_tolerance of {gap_tolerance}"
                        )
                        del temp_locus

    def prune_loci(self, min_genes=1):
        # remove loci if they don't have any features
        assembly_keys = list(self.assemblies.keys())
        for assembly_key in assembly_keys:
            loci_keys = list(self.assemblies[assembly_key].loci.keys())
            for loci_key in loci_keys:
                if (
                    len(self.assemblies[assembly_key].loci[loci_key].features)
                    < min_genes
                ):
                    del self.assemblies[assembly_key].loci[loci_key]

    def top_k_loci(self, topk=2):
        assembly_keys = list(self.assemblies.keys())
        for assembly_key in assembly_keys:
            loci_dict = {
                k: len(v.features)
                for k, v in self.assemblies[assembly_key].loci.items()
            }
            to_delete = [
                i
                for i in list(self.assemblies[assembly_key].loci.keys())
                if i not in sorted(loci_dict, key=loci_dict.get, reverse=True)[:topk]
            ]
            for loci_key in to_delete:
                del self.assemblies[assembly_key].loci[loci_key]

    ########################################
    # OUTPUTS FOR NEXTFLOW
    ########################################
    @staticmethod
    def _create_internal_locus_id(assembly_id, locus_id):
        # because locus id can't be assured to be unique across assemblies
        return hasher.sha512_hash(f"{assembly_id}___{locus_id}")

    def export_locus_to_protein(self, outdir: str = "."):
        """Locus to protein table for import into Neo4j

        Args:
            outdir (str, optional): Defaults to ".".
        """
        outpath = Path(outdir, "locus_to_protein")
        with open(outpath, "a") as handle:
            tsv_writer = csv.writer(
                handle, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL
            )
            for ak, av in self.assemblies.items():
                for k, loci in av.loci.items():
                    temp_list = list(loci.features)
                    # sort features by id then start to maintain consistent output
                    temp_list.sort(key=attrgetter("id"))
                    temp_list.sort(key=attrgetter("start"))
                    for feature in temp_list:
                        if feature.feature_is_protein():
                            # ["locus_id", "hash_id", "start", "end", "strand"]
                            tsv_writer.writerow(
                                [
                                    self._create_internal_locus_id(
                                        assembly_id=ak, locus_id=k
                                    ),
                                    feature.id,
                                    feature.start,
                                    feature.end,
                                    feature.strand,
                                ]
                            )

    def export_protein_info(self, outdir: str = "."):
        """Protein table for import into Neo4j

        Args:
            outdir (str, optional): Defaults to ".".
        """
        outpath = Path(outdir, "protein_info")
        with open(outpath, "a") as handle:
            tsv_writer = csv.writer(
                handle, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL
            )
            for protein in self.proteins.values():
                if protein.sequence is None:
                    prot_len = None
                else:
                    prot_len = len(protein.sequence)
                tsv_writer.writerow(
                    [
                        protein.hash_id,
                        # TODO: add database "source" of protein?
                        protein.other_id,
                        protein.description,
                        prot_len,  # protein length
                    ]
                )

    def export_assembly_to_locus(self, outdir: str = "."):
        """Assembly to locus table for import into Neo4j

        Args:
            outdir (str, optional): Defaults to ".".
        """
        outpath = Path(outdir, "assembly_to_locus")
        with open(outpath, "a") as handle:
            tsv_writer = csv.writer(
                handle, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL
            )
            for ak, av in self.assemblies.items():
                for k in av.loci.keys():
                    #  ["assembly", "internal_locus_id"]
                    tsv_writer.writerow(
                        [ak, self._create_internal_locus_id(assembly_id=ak, locus_id=k)]
                    )

    def export_loci(self, outdir: str = "."):
        """Loci table for import into Neo4j

        Args:
            outdir (str, optional): Defaults to ".".
        """
        outpath = Path(outdir, "loci")
        with open(outpath, "a") as handle:
            tsv_writer = csv.writer(
                handle, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL
            )
            for ak, av in self.assemblies.items():
                for k in av.loci.keys():
                    internal_id = self._create_internal_locus_id(
                        assembly_id=ak, locus_id=k
                    )
                    temp_v = []
                    for i in self.assemblies[ak].loci[k].info.values():
                        if isinstance(i, list):
                            if len(i) > 1:
                                temp_v.append(";".join(i))
                            else:
                                temp_v.append(i[0])
                        else:
                            temp_v.append(i)
                    out_list = [internal_id] + [k] + temp_v
                    #  ["internal_locus_id" "locus_id"]
                    tsv_writer.writerow(out_list)

    def export_assembly(self, outdir: str = "."):
        """Assembly table for import into Neo4j

        Args:
            outdir (str, optional): Defaults to ".".
        """
        outpath = Path(outdir, "assemblies")
        with open(outpath, "a") as handle:
            tsv_writer = csv.writer(
                handle, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL
            )
            for k in self.assemblies.keys():
                temp_v = []
                for i in self.assemblies[k].info.values():
                    if isinstance(i, list):
                        if len(i) > 1:
                            temp_v.append(";".join(i))
                        else:
                            temp_v.append(i[0])
                    else:
                        temp_v.append(i)
                tsv_writer.writerow([k] + temp_v)

    def export_assemby_to_taxid(self, outdir: str = "."):
        """Assembly to taxid table for import into Neo4j

        Args:
            outdir (str, optional): Defaults to ".".
        """
        outpath = Path(outdir, "assembly_to_taxid")
        with open(outpath, "a") as handle:
            tsv_writer = csv.writer(
                handle, delimiter="\t", quotechar='"', quoting=csv.QUOTE_MINIMAL
            )
            has_taxid = 0
            has_no_taxid = 0
            for k, v in self.assemblies.items():
                if v.taxid:
                    tsv_writer.writerow([k, v.taxid])
                    has_taxid += 1
                else:
                    has_no_taxid += 1
            log.info(
                f"{has_taxid} assemblies were connected to their taxon; {has_no_taxid} assemblies were not connected to their taxon; "
            )

    def drop_non_protein_features(self, return_removed=False):
        """Drop features that aren't proteins

        Args:
            return_removed (bool, optional): . Defaults to False.

        Returns:
            list: list of features that were removed (if return_removed=True)
        """
        report = {}
        for a_k, a_v in self.assemblies.items():
            report[a_k] = {}
            for l_k, l_v in a_v.loci.items():
                report[a_k][l_k] = {"retained": 0, "removed": 0}
                kept_features = set()
                for feature in l_v.features:
                    if feature.feature_is_protein():
                        # TODO: remove the need for deepcopy
                        kept_features.add(deepcopy(feature))
                        report[a_k][l_k]["retained"] += 1
                    else:
                        report[a_k][l_k]["removed"] += 1
                self.assemblies[a_k].loci[l_k].features = kept_features
        if return_removed:
            return report
        else:
            del report
