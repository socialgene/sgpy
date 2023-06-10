from pathlib import Path

import socialgene.utils.file_handling as fh
from socialgene.utils.logging import log


class Domtblout:
    def __init__(self):
        pass

    def parse_domtblout(self, input_path, hmmsearch_or_hmmscan="hmmsearch"):
        """Parse a HMMER domtblout file

        Args:
            input_path (str): path to file
            hmmsearch_or_hmmscan (str): Determines which acc belongs to protein and which belongs to the hmm

        Raises:
            FileNotFoundError: _description_
        """
        # For hmmsearch/hmmscan, HMMER switches the query/target for domtblout, so make that variable here
        for i in self._parse_domtblout(
            input_path=input_path, hmmsearch_or_hmmscan=hmmsearch_or_hmmscan
        ):
            _ = self.add_protein(hash_id=i["protein_id"])
            self.proteins[i["protein_id"]].add_domain(**i)

    def _parse_domtblout(self, input_path, hmmsearch_or_hmmscan="hmmsearch"):
        """Parse a HMMER domtblout file

        Args:
            input_path (str): path to file
            hmmsearch_or_hmmscan (str): Determines which acc belongs to protein and which belongs to the hmm

        Raises:
            FileNotFoundError: _description_
        """
        # For hmmsearch/hmmscan, HMMER switches the query/target for domtblout, so make that variable here
        if hmmsearch_or_hmmscan == "hmmscan":
            target_name = "hmm_id"
            query_name = "protein_id"
        elif hmmsearch_or_hmmscan == "hmmsearch":
            target_name = "protein_id"
            query_name = "hmm_id"
        else:
            raise ValueError
        all_columns = (
            target_name,
            "target_acc",
            "tlen",
            query_name,
            "query_acc",
            "qlen",
            "evalue",
            "seq_pro_score",
            "seq_pro_bias",
            "number",
            "of",
            "c_evalue",
            "i_evalue",
            "domain_score",
            "domain_bias",
            "hmm_from",
            "hmm_to",
            "ali_from",
            "ali_to",
            "env_from",
            "env_to",
            "acc",
            "desc",
        )

        input_path = Path(input_path)
        if not input_path.exists():
            raise FileNotFoundError(input_path)
        with fh.open_file(input_path) as f:
            for line in f:
                if line.startswith("#"):
                    pass
                else:
                    yield dict(zip(all_columns, line.split()))


class ParsedDomtblout:
    """For HMMER domtblout files already processed and rewritten with class Domtblout()"""

    def __init__(self):
        pass

    def parse_parseddomtblout(self, input_path):
        input_path = Path(input_path)
        with fh.open_file(input_path) as f:
            for line in f:
                line = line.strip().split("\t")
                self.add_protein(
                    hash_id=line[0],
                )
                self.proteins[line[0]].add_domain(
                    hash_id=line[1],
                    env_from=int(line[2]),
                    env_to=int(line[3]),
                    seq_pro_score=float(line[4]),  # maybe 7? TODO:
                    evalue=float(line[5]),
                    i_evalue=float(line[6]),
                    domain_bias=float(line[7]),  # or 8
                    domain_score=float(line[8]),  # or 7
                    seq_pro_bias=float(line[9]),  # ?
                    hmm_from=int(line[10]),
                    hmm_to=int(line[11]),
                    ali_from=int(line[12]),
                    ali_to=int(line[13]),
                )
        log.info(
            f"Read {str(len(self.proteins))} proteins and {str(sum([len(i.domains) for i in self.proteins.values()]))} domains from {str(input_path)}"
        )


class HmmerParser(Domtblout, ParsedDomtblout):
    def __init__(self):
        super().__init__()

    def parse_hmmout(self, filepath, **kwargs):
        """Parse sequence files (main function)
        Args:
            filepath (str): path to sequence file
        """
        # TODO: fh.check_if_tar(filepath=filepath)

        if check_if_domtblout(filepath):
            self.parse_domtblout(filepath, **kwargs)
        elif check_if_parseddomtblout(filepath):
            self.parse_parseddomtblout(**kwargs)
        else:
            raise NotImplementedError(
                "May not be implemented, or you need to use the genbank/fasta parser directly. (e.g. for tar archives)"
            )


def check_if_parseddomtblout(filepath):
    """Check if headers are maybe from a parsed domtblout file

    Args:
        filepath (str): path to sequence file

    Returns:
        bool: true/false
    """
    log.info(filepath)
    with fh.open_file(filepath) as f:
        l1 = f.readline()
        l2 = f.readline()
        a = l1.count("\t")
        b = l2.count("\t")
    return all(i == 13 for i in [a, b])


def check_if_domtblout(filepath):
    """Check if headers are from a HMMER domtblout file
    Args:
        filepath: file path of file to guess
    Returns:
        bool: true/false
    """
    with fh.open_file(filepath) as f:
        l1 = f.readline()
        l2 = f.readline()
    expected_header = "#---fullsequence-----------------thisdomain-------------hmmcoordalicoordenvcoord\n"
    if l1.replace(" ", "") == expected_header:
        if (
            l2.replace(" ", "")
            == "#targetnameaccessiontlenquerynameaccessionqlenE-valuescorebias#ofc-Evaluei-Evaluescorebiasfromtofromtofromtoaccdescriptionoftarget\n"
        ):
            return True
        else:
            return False
    else:
        return False
