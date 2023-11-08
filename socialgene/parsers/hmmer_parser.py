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
            self.add_protein(uid=i["external_id"])
            self.proteins[i["external_id"]].add_domain(**i)

    def _parse_domtblout(self, input_path, hmmsearch_or_hmmscan="hmmsearch"):
        """
        The function `_parse_domtblout` is used to parse a HMMER domtblout file, extracting relevant
        information and returning it as a dictionary.

        Args:
          input_path: The `input_path` parameter is a string that represents the path to the HMMER
        domtblout file that needs to be parsed.
          hmmsearch_or_hmmscan: The parameter `hmmsearch_or_hmmscan` is a string that determines whether
        the input file is from `hmmsearch` or `hmmscan`. Defaults to hmmsearch
        """
        # For hmmsearch/hmmscan, HMMER switches the query/target for domtblout, so make that variable here
        if hmmsearch_or_hmmscan == "hmmscan":
            target_name = "hmm_id"
            query_name = "external_id"
        elif hmmsearch_or_hmmscan == "hmmsearch":
            target_name = "external_id"
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
        with fh.open_read(input_path) as f:
            for line in f:
                if line.startswith("#"):
                    pass
                else:
                    yield dict(zip(all_columns, line.split()))


class ParsedDomtblout:
    """For HMMER domtblout files already processed and rewritten with class Domtblout()"""

    def __init__(self):
        pass

    def parse_parseddomtblout(self, input_path, **kwargs):
        """
        The function `parse_parseddomtblout` reads data from a parsed domtblout TSV and adds protein and domain
        information to a data structure.

        Args:
          input_path: The `input_path` parameter is the path to the input file that contains the data to
        be parsed. It should be a string representing the file path.
        """
        input_path = Path(input_path)
        with fh.open_read(input_path) as f:
            for line in f:
                line = line.strip().split("\t")
                self.add_protein(
                    uid=line[0],
                )
                print(bool(line[14]))
                self.proteins[line[0]].add_domain(
                    uid=line[1],
                    env_from=int(line[2]),
                    env_to=int(line[3]),
                    seq_pro_score=float(line[4]),
                    evalue=float(line[5]),
                    i_evalue=float(line[6]),
                    domain_bias=float(line[7]),
                    domain_score=float(line[8]),
                    seq_pro_bias=float(line[9]),
                    hmm_from=int(line[10]),
                    hmm_to=int(line[11]),
                    ali_from=int(line[12]),
                    ali_to=int(line[13]),
                    exponentialized=bool(line[14]),
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
        if check_if_domtblout(filepath):
            self.parse_domtblout(filepath, **kwargs)
        elif check_if_parseddomtblout(filepath):
            self.parse_parseddomtblout(**kwargs)
        else:
            raise NotImplementedError(
                "May not be implemented, or you need to use the genbank/fasta parser directly. (e.g. for tar archives)"
            )


def check_if_parseddomtblout(filepath):
    """
    The function `check_if_parseddomtblout` checks if the headers in a file are possibly from a parsed
    domtblout file.

    Args:
      filepath: The `filepath` parameter is a string that represents the path to a SocialGene parsed domtblout file.

    Returns:
      a boolean value, either True or False.
    """
    log.info(filepath)
    with fh.open_read(filepath) as f:
        l1 = f.readline()
        l2 = f.readline()
        a = l1.count("\t")
        b = l2.count("\t")
    return all(i == 14 for i in [a, b])


def check_if_domtblout(filepath):
    """
    The function `check_if_domtblout` checks if the headers of a file match the expected headers of a
    HMMER domtblout file.

    Args:
      filepath: The `filepath` parameter is a string that represents the file path of the file that
    needs to be checked.

    Returns:
      a boolean value (True or False) indicating whether the headers in the file at the given filepath
    match the expected headers for a HMMER domtblout file.
    """
    with fh.open_read(filepath) as f:
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
