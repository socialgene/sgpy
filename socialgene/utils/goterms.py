import argparse
import csv
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


def write_nodes(G, outpath):
    with open(outpath, "w") as h:
        writer = csv.writer(
            h,
            quoting=csv.QUOTE_MINIMAL,
            delimiter="\t",
        )
        for node in G.nodes:
            writer.writerow(
                (
                    node,
                    G.nodes[node].get("namespace", None),
                    G.nodes[node].get("name", None),
                    G.nodes[node].get("def", None),
                )
            )


def write_edges(G, outpath):
    with open(outpath, "w") as h:
        for edge in G.edges:
            h.write(f"{edge[0]}\t{edge[1]}\tGO_{edge[2].upper()}\n")


def main():
    args = parser.parse_args()
    log.info(f"Downloading and parsing {OBO_URL}")
    graph = download(url=OBO_URL)
    graph = nx.relabel_nodes(graph, lambda x: x.removeprefix("GO:"))
    log.info("Writing the nodes")
    write_nodes(graph, Path(args.outdir, "goterms"))
    log.info("Writing the edgelist")
    # data=False makes write_edgelist only return node pairs
    write_edges(graph, Path(args.outdir, "goterm_edgelist"))
    log.info(f"Finshed reading/writing {OBO_URL}")


if __name__ == "__main__":
    main()
