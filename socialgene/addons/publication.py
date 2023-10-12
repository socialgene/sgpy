import re


class Publication:
    __slots__ = [
        "doi",
        "pmid",
        "authors",
        "title",
        "journal",
        "year",
    ]

    def __init__(
        self,
        doi: str = None,
        pmid: str = None,
        authors: str = None,
        title: str = None,
        journal: str = None,
        year: str = None,
        **kwargs,
    ) -> None:
        self.doi = doi
        self.pmid = pmid
        self.authors = authors
        self.title = title
        self.journal = journal
        self.year = year

    @staticmethod
    def _extract_doi(input):
        return re.search("10.\\d{4,9}/[-._;()/:a-z0-9A-Z]+", input).group()
