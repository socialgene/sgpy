from collections import OrderedDict
from uuid import uuid4

import socialgene.hashing.hashing as hasher
from socialgene.config import env_vars
from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log
from socialgene.utils.simple_math import find_exp


class LocusAssemblyMetadata:
    # fmt: off
    __slots__= ["altitude","bio_material","bioproject","biosample","cell_line","cell_type","chromosome","clone","clone_lib","collected_by","collection_date","country","cultivar","culture_collection","db_xref","dev_stage","ecotype","environmental_sample","focus","germline","haplogroup","haplotype","host","identified_by","isolate","isolation_source","lab_host","lat_lon","macronuclear","map","mating_type","metagenome_source","mol_type","note","organelle","organism","pcr_primers","plasmid","pop_variant","proviral","rearranged","segment","serotype","serovar","sex","specimen_voucher","strain","sub_clone","submitter_seqid","sub_species","sub_strain","tissue_lib","tissue_type","transgenic","type_material","variety"]  # noqa: E231,E225
    # fmt: on

    def __init__(self, **kwargs) -> None:
        self.update(kwargs)

    @property
    def all_attributes(self):
        return OrderedDict({i: self.__getitem__(i) for i in self.__slots__})

    @property
    def values(self):
        return self.all_attributes.values()

    def __getitem__(self, item):
        try:
            return self.__getattribute__(item)
        except Exception as e:
            log.debug(e)
            return None

    # add a function that takes a dict as input and sets any attributes in the dict keys
    def update(self, d: dict):
        for k, v in d.items():
            temp_v = None
            if isinstance(v, list):
                # if info is a list, collapse into a single string
                if len(v) > 1:
                    temp_v = ";".join(v)
                else:
                    temp_v = v[0]
            else:
                temp_v = v
            temp_v.replace("\t", "")
            temp_v.replace("\n", "")
            temp_v.replace("\r", "")
            temp_v.strip()
            if k in self.__slots__:
                setattr(self, k, temp_v)


class ProteinSequence:
    """Class used for working with protein sequences and can be initialized with either a sequence or a hash_id."""

    _amino_acids = [
        "A",
        "R",
        "N",
        "D",
        "C",
        "Q",
        "E",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
        "X",
        "Z",
        "J",
        "U",
        "B",
        "O",
        "*",
    ]
    __slots__ = ["hash_id", "__crc64", "__md5", "sequence"]

    def __init__(self, sequence: str = None, hash_id: str = None):
        """Class for holding an amino acid sequence (protein)

        Args:
            sequence (str, optional): amino acid sequence
            hash_id (str, optional): a hash_id can be provided if a sequence isn't

        Raises:
            ValueError: Must provide either a sequence or a hash_id
        """
        # the sequence input variable has a default so that domtblout can create a socialgene object
        self.sequence = sequence
        if isinstance(self.sequence, type(None)):
            if isinstance(hash_id, type(None)):
                raise ValueError("Must provide either a sequence or a hash_id")
            else:
                self.hash_id = hash_id
        else:
            self._assign_hash()

    @property
    def all_attributes(self):
        """
        The function returns a dictionary containing the attributes of an object.

        Returns:
          The `__dict__` method is returning a dictionary that contains the names and values of all the
        attributes of the object. The attributes are obtained using the `getattr` function and are filtered
        to only include attributes that are defined in the `__slots__` list and that exist in the object.
        The dictionary is sorted based on the names of the attributes.
        """
        return {s: getattr(self, s) for s in sorted(self.__slots__) if hasattr(self, s)}

    def _one_letter_amino_acids(self):
        """Create an ordered dictionary of amino acids. Used to count AAs in a protein sequence.

        Returns:
            dict: amino acid count
        """
        return OrderedDict({i: 0 for i in self._amino_acids})

    def _amino_acid_count(self):
        """Create a '-' separated string of an amino acid count

        Returns:
            str: sequence of amino acids
        """
        if self.sequence is None:
            return "-".join(["0" for i in self._one_letter_amino_acids()])
        else:
            return "-".join(
                [str(self.sequence.count(i)) for i in self._one_letter_amino_acids()]
            )

    def _assign_hash(self):
        """
        The function assigns a hash value to a sequence of amino acids.
        """
        self._standardize_sequence()
        self.hash_id = hasher.hash_aminos(self.sequence)

    @property
    def crc64(self):
        return hasher.hash_aminos(self.sequence, algo="crc64")

    @property
    def md5(self):
        return hasher.hash_aminos(self.sequence, algo="md5")

    def _standardize_sequence(self):
        """
        The function converts the protein sequence to uppercase and checks if all characters are valid
        amino acids, raising an error if an unknown character is found.
        """
        self.sequence = self.sequence.upper()
        if not all([i in self._amino_acids for i in set(self.sequence)]):
            log.error(self.sequence)
            raise ValueError("Unknown character/letter in protein sequence")

    def sequence_length(self):
        """Return the length of the protein

        Returns:
            int: number of amino acids
        """
        return len(self.sequence)


class Location:
    __slots__ = ["start", "end", "strand"]
    # TODO: handle zero indexing

    def __init__(
        self,
        start: int = None,
        end: int = None,
        strand: int = None,
        **kwargs,  # this kwarg isn't accessed but is here so that calling Location with dict unpacking with extra args doesn't fail
    ):
        """Class describing genomic coordinates and strand direction.

        Args:
            start (int, optional): start coordinate
            end (int, optional): end coordinate
            strand (int, optional): DNA strand
        """
        self.start = start
        self.end = end
        self.strand = strand

    @property
    def all_attributes(self):
        return {s: getattr(self, s) for s in sorted(self.__slots__) if hasattr(self, s)}


class Domain:
    """Class for holding information about a domain/motif annotation"""

    __slots__ = [
        "hmm_id",
        "env_from",  # page 38 of the HMMER User Guide (http://eddylab.org/software/hmmer/Userguide.pdf) suggests using envelope (not hmm_from/to, ali_from/to)
        "env_to",
        "seq_pro_score",
        "evalue",
        "i_evalue",
        "domain_bias",
        "domain_score",
        "seq_pro_bias",
        "hmm_from",
        "hmm_to",
        "ali_from",
        "ali_to",
        "exponentialized",
    ]

    def __init__(
        self,
        hmm_id: str = None,
        env_from: int = None,
        env_to: int = None,
        seq_pro_score: float = None,
        evalue: float = None,
        i_evalue: float = None,
        domain_bias: float = None,
        domain_score: float = None,
        seq_pro_bias: float = None,
        hmm_from: int = None,
        hmm_to: int = None,
        ali_from: int = None,
        ali_to: int = None,
        exponentialized: bool = True,
        **kwargs,  # this kwarg isn't accessed but is here so that calling Domain with dict unpacking with extra args doesn't fail
    ):
        """
        Class for holding information about a domain/motif annotation

        Args:
          hmm_id (str): The `hmm_id` parameter is a string that represents the identifier of the hidden
        Markov model (HMM) associated with the domain.
          env_from (int): The `env_from` parameter represents the starting position of the domain in the
        target sequence.
          env_to (int): The `env_to` parameter represents the end position of the environment (sequence)
        in the domain. It is an integer value.
          seq_pro_score (float): The `seq_pro_score` parameter is a floating-point number that
        represents the sequence profile score of the domain.
          evalue (float): The `evalue` parameter is a floating-point number that represents the E-value
        of the domain. The E-value is a statistical measure that indicates the expected number of
        domains with a similar score or better that would occur by chance in a database of the same
        size.
          i_evalue (float): The `i_evalue` parameter is a floating-point number that represents the
        independent E-value of a domain. It is used to assess the statistical significance of the match
        between the domain and the sequence.
          domain_bias (float): The `domain_bias` parameter is a float that represents the bias score of
        a domain. It is used to measure the likelihood that a domain is present in a sequence due to
        chance rather than functional significance.
          domain_score (float): The `domain_score` parameter is a float that represents the score of the
        domain.
          seq_pro_bias (float): The `seq_pro_bias` parameter is a float value that represents the
        sequence profile bias of a domain. It is used to measure the bias in the sequence profile
        alignment.
          hmm_from (int): The `hmm_from` parameter represents the starting position of the domain in the
        HMM (Hidden Markov Model) profile. It is an integer value.
          hmm_to (int): The `hmm_to` parameter is an integer that represents the ending position of the
        hidden Markov model (HMM) alignment in the domain.
          ali_from (int): The `ali_from` parameter represents the starting position of the alignment in
        the sequence alignment. It is an integer value.
          ali_to (int): The `ali_to` parameter represents the ending position of the alignment in the
        sequence.
          exponentialized (bool): A boolean flag indicating whether the evalue and i_evalue should be
        exponentialized or not. If set to True, the evalue and i_evalue will be converted to exponential
        form. Defaults to True
        """
        super(Domain, self).__init__()
        if not isinstance(exponentialized, bool):
            raise ValueError(
                f"exponentialized mus be bool, was {type(exponentialized)}"
            )
        self.exponentialized = exponentialized
        if exponentialized:
            self.evalue = find_exp(evalue)
            self.i_evalue = find_exp(i_evalue)
        else:
            self.evalue = round(float(evalue), 1)
            self.i_evalue = round(float(i_evalue), 1)
        # some of are rounded because of differences between pyhmmer and hmmer results
        self.hmm_id = str(hmm_id)
        self.env_from = int(env_from)
        self.env_to = int(env_to)
        self.seq_pro_score = round(float(seq_pro_score), 1)
        self.domain_bias = round(float(domain_bias), 1)
        self.domain_score = round(float(domain_score), 1)
        self.seq_pro_bias = round(float(seq_pro_bias), 1)
        self.hmm_from = int(hmm_from)
        self.hmm_to = int(hmm_to)
        self.ali_from = int(ali_from)
        self.ali_to = int(ali_to)

    @property
    def all_attributes(self):
        ord_dict = OrderedDict()
        ord_dict["hmm_id"] = str(self.hmm_id)
        ord_dict["env_from"] = int(self.env_from)
        ord_dict["env_to"] = int(self.env_to)
        ord_dict["seq_pro_score"] = round(self.seq_pro_score, 1)
        ord_dict["evalue"] = round(self.evalue, 1)
        ord_dict["i_evalue"] = round(self.i_evalue, 1)
        ord_dict["domain_bias"] = round(self.domain_bias, 1)
        ord_dict["domain_score"] = round(self.domain_score, 1)
        ord_dict["seq_pro_bias"] = round(self.seq_pro_bias, 1)
        ord_dict["hmm_from"] = int(self.hmm_from)
        ord_dict["hmm_to"] = int(self.hmm_to)
        ord_dict["ali_from"] = int(self.ali_from)
        ord_dict["ali_to"] = int(self.ali_to)
        ord_dict["exponentialized"] = bool(self.exponentialized)
        return ord_dict

    def get_hmm_id(self):
        return self.hmm_id

    @property
    def domain_within_threshold(self):
        """Check if a protein's domain annotation is within the threshold for inclusion (currently i_evalue based)
        Returns:
            bool: is less than or equal to the set ievalue cutoff?
        """
        if self.exponentialized:
            return self.i_evalue <= find_exp(env_vars["HMMSEARCH_IEVALUE"])
        else:
            return self.i_evalue <= env_vars["HMMSEARCH_IEVALUE"]

    def __hash__(self):
        """Used to prevent adding duplicate domains
        Returns:
            hash: hash for set()
        """
        return hash((self.hmm_id, self.env_from, self.env_to, self.i_evalue))

    def __eq__(self, other):  # pragma: no cover
        if not isinstance(other, type(self)):
            return NotImplemented
        return (
            self.hmm_id == other.hmm_id
            and self.env_from == other.env_from
            and self.env_to == other.env_to
            and self.i_evalue == other.i_evalue
        )


class Protein(
    ProteinSequence,
):
    """Container class for describing a single protein"""

    # seqlen is included as a variable so it can be provided (e.g. a sequence isn't provided and proteins are created manually)
    __slots__ = ["description", "external_protein_id", "seqlen", "domains"]

    def __init__(
        self,
        description: str = None,
        external_protein_id: str = None,
        seqlen: int = None,
        domains: set = None,
        *args,
        **kwargs,
    ):
        """_summary_

        Args:
            description (str, optional): Protein description.
            external_protein_id (str, optional): Non-hash-id descriptor (usually a database accession, e.g. NCBI's).
            seqlen (int, optional): Amino acid sequence length
            domains (Set, optional): Set of Domain() objects.
        """
        super().__init__(*args, **kwargs)
        self.description = description
        self.external_protein_id = external_protein_id
        self.seqlen = seqlen
        self.domains = domains if domains is not None else set()

    @property
    def all_attributes(self):
        return {s: getattr(self, s) for s in sorted(self.__slots__) if hasattr(self, s)}

    def add_domain(
        self,
        *args,
        **kwargs,
    ):
        """Add a domain to a protein if it meets threshold (ievalue)

        Returns:
            bool: used for counting # of domains added to protein
        """
        self.domains.add(Domain(*args, **kwargs))

    def add_domains_from_neo4j(self):
        try:
            with GraphDriver() as db:
                for i in db.run(
                    """
                    MATCH (p1:protein)<-[a1:ANNOTATES]-(h1:hmm)
                    where p1.uid in ["_xvHFg1H1f43WxPcZT1P8Qbcx60chSxh"]
                    RETURN apoc.map.merge({hmm_id:h1.uid}, properties(a1))
                    """,
                    uid=str(self.hash_id),
                ):
                    self.add_domain(**i.value())
        except Exception:
            log.debug(f"Error trying to retrieve domains for {self.hash_id}")

    @property
    def domain_list_sorted_by_mean_envelope_position(self):
        # (ie can't sort a set())
        # 'if' is to check whether domains is None
        # A tuple is passed so that the HMM annnotations are first sorted by the envelope midpoint
        # then ties are sorted alphabetically based on hmm_id for consistency
        return sorted(
            list(self.domains),
            key=lambda tx: ((tx.env_from + tx.env_to) / 2, tx.hmm_id),
        )

    @property
    def domain_vector(
        self,
    ):
        """Get the domain hash_ids for a protein as an ordered list

        list: list of domain hash_ids
        """
        if not self.domains:
            log.debug(f"Tried to get domains from domain-less protein {self}")
            return []
        return [
            i.get_hmm_id() for i in self.domain_list_sorted_by_mean_envelope_position
        ]

    def filter_domains(self):
        """Prune all domains in all proteins that don't meet the inclusion threshold (currently HMMER's i_evalue)"""
        _before_count = len(self.domains)
        temp = []
        self.domains = {i for i in self.domains if i.domain_within_threshold}

        del temp
        log.debug(
            f"Removed {str(_before_count - len(self.domains))} domains from {self.external_protein_id}"
        )

    @property
    def fasta_string_defline_hash_id(self):
        return f">{self.hash_id}\n{self.sequence}\n"

    @property
    def fasta_string_defline_external_id(self):
        return f">{self.external_protein_id}\n{self.sequence}\n"


class Feature(Location):
    """Container class for describing a feature on a locus"""

    __slots__ = [
        "parent_object",
        "protein_hash",
        "protein_id",
        "type",
        "locus_tag",
        "description",
        "note",
        "goterms",
        "partial_on_complete_genome",
        "missing_start",
        "missing_stop",
        "internal_stop",
        "partial_in_the_middle_of_a_contig",
        "missing_N_terminus",
        "missing_C_terminus",
        "frameshifted",
        "too_short_partial_abutting_assembly_gap",
        "incomplete",
    ]

    def __init__(
        self,
        parent_object=None,
        protein_hash: str = None,
        protein_id: str = None,
        type: str = None,
        locus_tag: str = None,
        description=None,
        note=None,
        goterms=None,
        partial_on_complete_genome=None,
        missing_start=None,
        missing_stop=None,
        internal_stop=None,
        partial_in_the_middle_of_a_contig=None,
        missing_N_terminus=None,
        missing_C_terminus=None,
        frameshifted=None,
        too_short_partial_abutting_assembly_gap=None,
        incomplete=None,
        **kwargs,
    ):
        """
        The above function is a constructor for a class that represents a feature on a locus, with
        various attributes and optional arguments.

        Args:
          protein_hash (str): A string representing the hash value of the protein.
          protein_id (str): The unique identifier for the protein associated with the feature.
          type (str): The "type" parameter is used to specify the type of feature. In this case, it is
        used to specify the type of the feature on a locus, such as "protein".
          locus_tag (str): The locus_tag parameter is a string that represents the unique identifier for
        a feature on a locus. It is typically used in genomics to identify a specific gene or protein
        within a genome.
          description: A description of the feature on a locus.
          note: The `note` parameter is used to provide additional information or comments about the
        feature. It is an optional parameter and can be any string value.
          goterms: The `goterms` parameter is used to store the Gene Ontology terms associated with the
        feature. Gene Ontology (GO) terms are standardized terms used to describe the function,
        location, and involvement of genes and gene products in biological processes.
          partial_on_complete_genome: A boolean flag indicating whether the feature is partial on a
        complete genome.
          missing_start: A boolean indicating whether the start of the feature is missing.
          missing_stop: A boolean flag indicating whether the stop codon of the feature is missing.
          internal_stop: A boolean flag indicating whether the feature has an internal stop codon.
          partial_in_the_middle_of_a_contig: This parameter is used to indicate whether the feature is
        partial and located in the middle of a contig.
          missing_N_terminus: This parameter is used to indicate whether the N-terminus (the starting
        end) of the feature is missing. It is a boolean value, where True indicates that the N-terminus
        is missing and False indicates that it is not missing.
          missing_C_terminus: The parameter "missing_C_terminus" is used to indicate whether the feature
        is missing the C-terminus (the end) of the protein sequence. It is a boolean value that can be
        set to True or False.
          frameshifted: The "frameshifted" parameter is a boolean flag that indicates whether the
        feature has a frameshift mutation. A frameshift mutation occurs when the reading frame of a gene
        is disrupted by the insertion or deletion of nucleotides, causing a shift in the codon reading
        frame and potentially altering the amino
          too_short_partial_abutting_assembly_gap: This parameter is used to indicate whether the
        feature is too short and abuts an assembly gap.
          incomplete: A boolean flag indicating whether the feature is incomplete.
        """
        super().__init__(**kwargs)
        self.parent_object = parent_object
        self.protein_hash = protein_hash
        self.protein_id = protein_id
        self.type = type
        self.locus_tag = locus_tag
        self.description = description
        self.note = note
        self.goterms = goterms
        self.partial_on_complete_genome = partial_on_complete_genome
        self.missing_start = missing_start
        self.missing_stop = missing_stop
        self.internal_stop = internal_stop
        self.partial_in_the_middle_of_a_contig = partial_in_the_middle_of_a_contig
        self.missing_N_terminus = missing_N_terminus
        self.missing_C_terminus = missing_C_terminus
        self.frameshifted = frameshifted
        self.too_short_partial_abutting_assembly_gap = (
            too_short_partial_abutting_assembly_gap
        )
        self.incomplete = incomplete

    @property
    def all_attributes(self):
        return {s: getattr(self, s) for s in sorted(self.__slots__) if hasattr(self, s)}

    def feature_is_protein(self):
        """Check if the feature "is a protein"

        Returns:
            bool: True if socialgene thinks it's a protein, False if not
        """
        # types of genbank file features to classify as a "protein"
        types_list = ["protein", "CDS"]
        return any([True for i in types_list if i == self.type])

    def __hash__(self):
        """Used to prevent adding duplicate features to a locus (for hash in set() in Assembly.add_locus())

        Returns:
            hash: hash for set()
        """
        return hash((self.end, self.start, self.protein_hash, self.strand, self.type))

    def __eq__(self, other):
        """Used for set() in Assembly.add_locus()"""
        if not isinstance(other, type(self)):
            return NotImplemented
        return (
            self.end == other.end
            and self.start == other.start
            and self.protein_hash == other.protein_hash
            and self.strand == other.strand
            and self.type == other.type
        )


class Locus:
    """Container holding a set() of genomic features"""

    __slots__ = ["parent_object", "features", "metadata", "external_id", "uid"]

    def __init__(self, parent_object=None, external_id=None):
        super().__init__()
        self.parent_object = parent_object
        self.features = set()
        self.metadata = LocusAssemblyMetadata()
        self.external_id = external_id
        self.uid = self.calc_uid()

    def add_feature(self, **kwargs):
        """Add a feature to a locus"""
        self.features.add(Feature(parent_object=self, **kwargs))

    def calc_uid(self):
        return hasher.hasher(f"{self.parent_object.uid}___{self.external_id}")

    @property
    def all_attributes(self):
        return {s: getattr(self, s) for s in sorted(self.__slots__) if hasattr(self, s)}

    @property
    def protein_hash_set(self):
        return {i.protein_hash for i in self.features}

    @property
    def features_sorted_by_midpoint(self):
        """Sorts features by mid-coordinate of each feature"""
        for i in sorted(list(self.features), key=lambda i: int((i.end + i.start) / 2)):
            yield i

    def drop_non_protein_features(self):
        self.features = {i for i in self.features if i.type == "protein"}


class Taxonomy:
    "Class is a reserved word so just underscore all ranks to be consistent"
    __slots__ = [
        "species_",
        "genus_",
        "family_",
        "order_",
        "class_",
        "phylum_",
        "clade_",
        "superkingdom_",
    ]

    def __init__(
        self,
        species_: str = None,
        genus_: str = None,
        family_: str = None,
        order_: str = None,
        class_: str = None,
        phylum_: str = None,
        clade_: str = None,
        superkingdom_: str = None,
        **kwargs,
    ) -> None:
        self.species_ = species_
        self.genus_ = genus_
        self.family_ = family_
        self.order_ = order_
        self.class_ = class_
        self.phylum_ = phylum_
        self.clade_ = clade_
        self.superkingdom_ = superkingdom_


class Assembly:
    """Container class holding a dictionary of loci (ie genes/proteins)"""

    __slots__ = ["loci", "taxid", "metadata", "uid", "taxonomy", "name"]

    def __init__(self, uid):
        super().__init__()
        self.uid = uid
        self.loci = {}
        self.taxid = None
        self.taxonomy = Taxonomy()
        self.metadata = LocusAssemblyMetadata()
        self.name = uid

    @property
    def protein_hash_set(self):
        """Return all protein hashes within assembly
        Returns:
            set: protein hashes
        """
        return set().union(*[i.protein_hash_set for i in self.loci.values()])

    @property
    def all_attributes(self):
        return {s: getattr(self, s) for s in sorted(self.__slots__) if hasattr(self, s)}

    def add_locus(self, external_id: str = None):
        """Add a locus to an assembly object"""
        if external_id is None:
            external_id = str(uuid4())
        if external_id not in self.loci:
            self.loci[external_id] = Locus(parent_object=self, external_id=external_id)
        else:
            log.debug(f"{external_id} already present")

    def fill_taxonomy_from_db(self):
        try:
            with GraphDriver() as db:
                res = db.run(
                    """
                        MATCH (a1:assembly {uid: $uid})-[:IS_TAXON]->(t2:taxid)-[:TAXON_PARENT*1..]->(t1:taxid)
                        WHERE t1.rank in ["genus", "family", "order", "class", "phylum", "clade", "superkingdom"]
                        WITH  apoc.map.fromLists([t2.rank],[t2.name]) as a, apoc.map.fromLists([t1.rank],[t1.name]) as b
                        return  apoc.map.mergeList(collect(a)+collect(b)) as tax_dict
                        """,
                    uid=str(self.uid),
                ).value()
                # Dict comprehension appends '_' to the keys, to match the args in Taxonomy
                self.taxonomy = Taxonomy(**{f"{k}_": v for k, v in res[0].items()})
        except Exception:
            log.debug(f"Error trying to retrieve taxonomy for {self.uid}")

    def fill_properties(self):
        try:
            with GraphDriver() as db:
                res = db.run(
                    """
                        MATCH (a1:assembly {uid: $uid})
                        return properties(a1)
                        """,
                    uid=self.uid,
                ).value()[0]
                self.metadata.update(res)
        except Exception:
            log.debug(f"Error trying to retrieve taxonomy for {self.uid}")


class Molbio:
    """Class for inheriting by SocialGene()"""

    def __init__(self):
        super().__init__()  # TODO: ?

        self.assemblies = {}  # TODO: ?
        self.proteins = {}  # TODO: ?

    def get_all_protein_hashes(self):
        """Return a list of all proteins hash_ids

        Returns:
            list: List of all proteins hash_ids
        """
        return [i.hash_id for i in self.proteins.values()]

    def add_protein(
        self,
        return_uid=True,
        **kwargs,
    ):
        """Add a protein to the protein dictionary

        Args:
            no_return (bool, optional): Whether the function should return the protein's hash_id. Defaults to False.

        Returns:
            str: Protein's hash
        """
        temp_protein = Protein(**kwargs)
        # only add protein if it doesn't already exist
        if temp_protein.hash_id not in self.proteins:
            # deepcopy teo ensure instances aren't shared
            self.proteins[temp_protein.hash_id] = temp_protein
        if return_uid:
            return temp_protein.hash_id

    def add_assembly(self, uid: str = None):
        """Add an assembly to a SocialGene object

        Args:
            uid (str, optional): Assembly identifier, should be unique across parsed input. Defaults to None.
        """
        if uid is None:
            uid = str(uuid4())
        if uid not in self.assemblies:
            self.assemblies[uid] = Assembly(uid)
        else:
            log.debug(f"{uid} already present")
