# python dependencies
from uuid import uuid4
from collections import OrderedDict

# external dependencies

# internal dependencies
from socialgene.config import env_vars
import socialgene.hashing.hashing as hasher
from socialgene.utils.logging import log
from socialgene.utils.simple_math import find_exp


class ProteinSequence:
    """Class for working with protein sequences"""

    amino_acids = [
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
    __slots__ = ["hash_id", "sequence"]

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
            self.hash()

    def one_letter_amino_acids(self):
        """Create an ordered dictionary of amino acids. Used to count AAs in a protein sequence.

        Returns:
            dict: amino acid count
        """
        return OrderedDict({i: 0 for i in self.amino_acids})

    def count_amino_acids(self):
        """Create a '-' separated string of an amino acid count

        Returns:
            str: sequence of amino acids
        """
        if self.sequence is None:
            return "NA"
        else:
            return "-".join(
                [str(self.sequence.count(i)) for i in self.one_letter_amino_acids()]
            )

    def hash(self):
        """Hash the amino acid sequence"""
        self._standardize_sequence()
        self.hash_id = hasher.hash_aminos(self.sequence)

    def _standardize_sequence(self):
        """Check an input AA sequence for expected characters

        Raises:
            ValueError: "Unknown character/letter in protein sequence
        """
        self.sequence = self.sequence.upper()
        if not all([i in self.amino_acids for i in set(self.sequence)]):
            log.error(self.sequence)
            raise ValueError("Unknown character/letter in protein sequence")

    def sequence_length(self):
        """Return the length of the protein

        Returns:
            int: number of amino acids
        """
        return len(self.sequence)


class Location:
    """Class describing genomic coordinates"""

    __slots__ = ["start", "end", "strand"]
    # TODO: handle zero indexing

    def __init__(
        self,
        start: int = None,
        end: int = None,
        strand: int = None,
    ):
        """Class describing genomic coordinates

        Args:
            start (int, optional): start coordinate
            end (int, optional): end coordinate
            strand (int, optional): DNA strand
        """
        self.start = start
        self.end = end
        self.strand = strand


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
    ]

    def __init__(
        self,
        exponentialized: bool = False,
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
        **kwargs,  # this kwarg isn't accessed but is here so that calling Domain with dict unpacking with extra args doesn't fail
    ):
        """Class for holding information about a domain/motif annotation

        Args:
            hmm_id (str, optional): hmm model hash id
            env_from (int, optional): see hmmer documentation for description
            env_to (int, optional): see hmmer documentation for description
            seq_pro_score (float, optional): see hmmer documentation for description
            e_value (float, optional): see hmmer documentation for description
            evalue (float, optional): see hmmer documentation for description
            domain_bias (float, optional): see hmmer documentation for description
            domain_score (float, optional): see hmmer documentation for description
            seq_pro_bias (float, optional): see hmmer documentation for description
            hmm_from (int, optional): see hmmer documentation for description
            hmm_to (int, optional): see hmmer documentation for description
            ali_from (int, optional): see hmmer documentation for description
            ali_to (int, optional): see hmmer documentation for description
        """
        super(Domain, self).__init__()
        # some of these are rounded because of differences between pyhmmer and hmmer results
        self.hmm_id = str(hmm_id)
        self.env_from = int(env_from)
        self.env_to = int(env_to)
        self.seq_pro_score = round(float(seq_pro_score), 1)
        if exponentialized:
            self.evalue = int(evalue)
            self.i_evalue = int(i_evalue)
        else:
            self.evalue = find_exp(evalue)
            self.i_evalue = find_exp(i_evalue)
        self.domain_bias = round(float(domain_bias), 1)
        self.domain_score = round(float(domain_score), 1)
        self.seq_pro_bias = round(float(seq_pro_bias), 1)
        self.hmm_from = int(hmm_from)
        self.hmm_to = int(hmm_to)
        self.ali_from = int(ali_from)
        self.ali_to = int(ali_to)

    def get_hmm_id(self):
        return self.hmm_id

    def domain_within_threshold(self):
        """Check if a protein's domain annotation is within the threshold for inclusion (currently i_evalue based)

        Raises:
            TypeError: Incorrect type for i_evalue

        Returns:
            int: exponent of the i_evalue
        """
        if not isinstance(self.i_evalue, int):
            raise TypeError(
                f"Domain i_evalue must be an integer representing the exponent, not {self.i_evalue}"
            )
        return self.i_evalue <= find_exp(float(env_vars["HMMSEARCH_IEVALUE"]))

    def get_dict(self):
        ord_dict = OrderedDict()
        ord_dict["hmm_id"] = str(self.hmm_id)
        ord_dict["env_from"] = int(self.env_from)
        ord_dict["env_to"] = int(self.env_to)
        ord_dict["seq_pro_score"] = round(self.seq_pro_score, 2)
        ord_dict["evalue"] = round(self.evalue, 2)
        ord_dict["i_evalue"] = round(self.i_evalue, 2)
        ord_dict["domain_bias"] = round(self.domain_bias, 2)
        ord_dict["domain_score"] = round(self.domain_score, 2)
        ord_dict["seq_pro_bias"] = round(self.seq_pro_bias, 2)
        ord_dict["hmm_from"] = int(self.hmm_from)
        ord_dict["hmm_to"] = int(self.hmm_to)
        ord_dict["ali_from"] = int(self.ali_from)
        ord_dict["ali_to"] = int(self.ali_to)
        return ord_dict

    def __hash__(self):
        """Used to prevent adding duplicate domains
        Returns:
            hash: hash for set()
        """
        return hash((self.hmm_id, self.env_from, self.env_to, self.i_evalue))

    def __eq__(self, other):
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

    __slots__ = ["description", "other_id", "seqlen", "domains"]

    def __init__(
        self,
        description: str = None,
        other_id: str = None,
        seqlen: int = None,
        domains: set = None,
        *args,
        **kwargs,
    ):
        """_summary_

        Args:
            description (str, optional): Protein description.
            other_id (str, optional): Non-hash-id descriptor (usually a database accession, e.g. NCBI's).
            seqlen (int, optional): Amino acid sequence length
            domains (Set, optional): Set of Domain() objects.
        """
        super().__init__(*args, **kwargs)
        self.description = description
        self.other_id = other_id
        self.seqlen = seqlen
        self.domains = domains if domains is not None else set()

    def only_id(self):
        self.description = None
        self.other_id = None
        self.domains = set()

    def add_domain(
        self,
        *args,
        **kwargs,
    ):
        """Add a domain to a protein if it meets threshold (ievalue)

        Returns:
            bool: used for counting # of domains added to protein
        """
        temp = Domain(*args, **kwargs)
        if temp.domain_within_threshold():
            self.domains.add(temp)
        # if temp.i_evalue is None:
        #     self.domains.append(deepcopy(temp))
        #     return True
        # if temp.domain_within_threshold():
        #     self.domains.append(deepcopy(temp))
        #     return True
        # else:
        #     return False

    def sort_domains_by_mean_envelope_position(self):
        # (ie can't sort a set())
        # 'if' is to check whether domains is None
        # A tuple is passed so that the HMM annnotations are first sorted by the envelope midpoint
        # then ties are sorted alphabetically based on hmm_id for consistency
        return sorted(
            list(self.domains),
            key=lambda tx: ((tx.env_from + tx.env_to) / 2, tx.hmm_id),
        )

    def get_domain_vector(
        self,
        only_unique=False,
    ):
        """Get the domain hash_ids for a protein as an ordered list

        Args:
            only_unique (bool, optional): Whether duplicate domain hash_ids should be included in the returned results. Defaults to False.

        Returns:
            list: list of domain hash_ids
        """
        if not self.domains:
            log.debug(f"Tried to get domains from domain-less protein {self}")
            return []
        if only_unique:
            # self.domains is a set() but that takes in to account more than just
            # the id so we need to get unique ids here
            return list(set([i.get_hmm_id() for i in self.domains]))
        else:
            sorted_domain_list = self.sort_domains_by_mean_envelope_position()
            return [i.get_hmm_id() for i in sorted_domain_list]


class Feature(Location):
    """Container class for describing a feature on a locus"""

    __slots__ = ["id", "type"]

    def __init__(self, id: str = None, type: str = None, **kwargs):
        """Container class for describing a feature on a locus

        Args:
            id (str, optional): id for the locus
            type (str, optional): type of feature e.g. "protein"
        """
        super().__init__(**kwargs)
        self.id = id
        self.type = type

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
        return hash((self.end, self.start, self.id, self.strand, self.type))

    def __eq__(self, other):
        """Used for set() in Assembly.add_locus()"""
        if not isinstance(other, type(self)):
            return NotImplemented
        return (
            self.end == other.end
            and self.start == other.start
            and self.id == other.id
            and self.strand == other.strand
            and self.type == other.type
        )


class Locus:
    """Container holding a set() of genomic features"""

    __slots__ = "features"

    def __init__(self):
        super().__init__()
        self.features = set()

    def add_feature(self, **kwargs):
        """Add a feature to a locus"""
        self.features.add(Feature(**kwargs))

    def sort_by_middle(self):
        """Sorts features by mid-coordinate of each feature"""
        self.features = list(self.features)
        self.features.sort(key=lambda i: int((i.end + i.start) / 2))


class Assembly:
    """Container class holding a dictionary of loci (ie genes/proteins)"""

    __slots__ = ["loci", "taxid"]

    def __init__(self):
        super().__init__()
        self.loci = {}
        self.taxid = None

    def add_locus(self, id: str = None):
        """Add a locus to an assembly object"""
        if id is None:
            id = str(uuid4())
        if id not in self.loci:
            self.loci[id] = Locus()

    def get_min_maxcoordinates(self):
        # TODO: broken
        """Get the minimum and maximum start/stop for a set of loci"""
        return {
            k: (min([i.start for i in v]), max([i.start for i in v]))
            for k, v in self.loci.items()
        }


class Molbio(Assembly):
    """Class for inheriting by SocialGene()"""

    def __init__(self):
        pass

        self.assemblies = {}
        self.proteins = {}

    def get_all_protein_hashes(self):
        """Return a list of all proteins hash_ids

        Returns:
            list: List of all proteins hash_ids
        """
        return [i.hash_id for i in self.proteins.values()]

    def add_protein(
        self,
        no_return=False,
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
        if not no_return:
            return temp_protein.hash_id

    def add_assembly(self, id: str = None, verbose=False):
        """Add an assembly to a SocialGene object

        Args:
            id (str, optional): Assembly identifier, should be unique across parsed input. Defaults to None.
            verbose (bool, optional): Print extra logging? Defaults to False.
        """
        if id is None:
            id = str(uuid4())
        if id not in self.assemblies:
            self.assemblies[id] = Assembly()
        elif verbose:
            log.debug(f"{id} already present")
