from socialgene.neo4j.neo4j import GraphDriver
from socialgene.utils.logging import log


class Substrate:
    # __slots__ = ["name", "proteinogenic", "smiles"]
    def __init__(self, name=None, proteinogenic=None, structure=None, **kwargs) -> None:
        self.name = name
        self.proteinogenic = proteinogenic
        self.smiles = structure

    def get(self):
        # only return dict key,values that aren't none
        # because won't be filtered in Neo4j step
        return {k: v for k, v in self.__dict__.items() if v}


class PKS:
    __slots__ = [
        "assembly",
        "gene_id",
        "cyclic",
        "release_type",
        "starter_unit",
        "subclasses",
        "synthases",
    ]

    def __init__(self) -> None:
        self.assembly = None
        self.gene_id = None
        self.cyclic = None
        self.release_type = None
        self.starter_unit = None
        self.subclasses = None
        self.synthases = None

    # def _assign_subclass():
    #     try:
    #         input_json["cluster"]["polyketide"]["synthases"][0]["subclass"]
    #     except Exception:
    #         pass

    def assign(self, input_json):
        self.assembly = input_json["cluster"]["mibig_accession"]
        if "nrp" in input_json["cluster"]:
            if "cyclic" in input_json["cluster"]["nrp"]:
                self.cyclic = input_json["cluster"]["nrp"]["cyclic"]
            if "nrps_genes" in input_json["cluster"]["nrp"]:
                for gene in input_json["cluster"]["nrp"]["nrps_genes"]:
                    self.gene_id = gene["gene_id"]
                    if "modules" in gene:
                        for module in gene["modules"]:
                            if (
                                "a_substr_spec" in module
                                and "substrates" in module["a_substr_spec"]
                            ):
                                self.module_type = "A_SUBSTR_SPEC"
                                if "publications" in list(module.values())[0]:
                                    self.publications = "; ".join(
                                        list(module.values())[0]["publications"]
                                    )
                                if "evidence" in list(module.values())[0]:
                                    self.evidence = "; ".join(
                                        list(module.values())[0]["evidence"]
                                    )
                                for substrate in list(module.values())[0]["substrates"]:
                                    self.substrates.append(Substrate(**substrate).get())

    def write_to_neo4j(self):
        if not self.gene_id:
            log.info("No NRPS")
        else:
            summary = (
                GraphDriver()
                .driver.execute_query(
                    f"""
                    MATCH ()-[e1:ENCODES {{protein_id:$gene_id}}]->(p1:protein)
                    UNWIND $substrates as substrate
                    MERGE (sub:substrate {{name: substrate.name}})
                    SET sub += substrate
                    MERGE (p1)-[r:{self.module_type} ]->(sub)
                    SET r += $relationship_info
                    """,
                    relationship_info={
                        "publications": self.publications,
                        "evidence": self.evidence,
                    },
                    substrates=self.substrates,
                    gene_id=self.gene_id,
                    database_="neo4j",
                )
                .summary
            )
            if summary.metadata.get("stats"):
                log.info(
                    f"{summary.metadata.get('stats').get('properties-set')} properties modified"
                )
            else:
                log.info("No properties modified")
        if self.cyclic:
            summary = (
                GraphDriver()
                # using "cyclic_nrps" becuse not sure if cyclic is specified
                # independently in mibig for pks and nrps
                .driver.execute_query(
                    """
                    MATCH (a1:assembly {$assembly})
                    SET a1.cyclic_nrps = $cyclic
                    """,
                    assembly=self.assembly,
                    cyclic=self.cyclic,
                    database_="neo4j",
                ).summary
            )
            if summary.metadata.get("stats"):
                log.info(
                    f"{summary.metadata.get('stats').get('properties-set')} properties modified"
                )
            else:
                log.info("No properties modified")
