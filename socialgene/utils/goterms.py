import argparse
from pathlib import Path

import networkx as nx
import obonet

from socialgene.utils.logging import log

parser = argparse.ArgumentParser(
    description="Download and create edgelist of all Go terms"
)

parser.add_argument(
    "--outdir",
    metavar="filepath",
    help="Output directory filepath",
    required=True,
)

OBO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"


def download(url):
    return obonet.read_obo(url)


def main():
    args = parser.parse_args()
    log.info(f"Downloading and parsing {OBO_URL}")
    graph = download(url=OBO_URL)
    graph = nx.relabel_nodes(graph, lambda x: x.removeprefix("GO:"))
    log.info("Writing the nodes")
    with open(Path(args.outdir, "goterms"), "w") as h:
        for i in sorted(set(graph.nodes.keys())):
            h.write(f"{i}\n")
    log.info("Writing the edgelist")
    # data=False makes write_edgelist only return node pairs
    nx.write_edgelist(
        graph, Path(args.outdir, "goterm_edgelist"), data=False, delimiter="\t"
    )
    log.info(f"Finshed reading/writing {OBO_URL}")


if __name__ == "__main__":
    main()
