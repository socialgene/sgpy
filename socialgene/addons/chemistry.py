from socialgene.base.chem import ChemicalCompound
from socialgene.neo4j.neo4j_element import Node, Relationship


class ChemicalFragment(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="chemical_fragment",
            description="Represents a chemical fragment as defined by rdkit.Chem.Descriptors",
            required_properties=["uid"],
            properties={
                "uid": str,
            },
        )


class ChemicalCompoundNode(Node):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="chemical_compound",
            description="Represents a chemical compound",
            required_properties=["inchi", "CanonicalSmiles"],
            properties={
                "uid": str,
                "MolWt": float,
                "HeavyAtomMolWt": float,
                "ExactMolWt": float,
                "NumValenceElectrons": int,
                "NumRadicalElectrons": int,
                "HeavyAtomCount": int,
                "NumAliphaticCarbocycles": int,
                "NumAliphaticHeterocycles": int,
                "NumAliphaticRings": int,
                "NumAromaticCarbocycles": int,
                "NumAromaticHeterocycles": int,
                "NumAromaticRings": int,
                "NumHAcceptors": int,
                "NumHDonors": int,
                "NumHeteroatoms": int,
                "NumRotatableBonds": int,
                "NumSaturatedCarbocycles": int,
                "NumSaturatedHeterocycles": int,
                "NumSaturatedRings": int,
                "RingCount": int,
                "MolLogP": float,
                "MolMR": float,
                "AnonymousGraph": str,
                "ElementGraph": str,
                "CanonicalSmiles": str,
                "MurckoScaffold": str,
                "ExtendedMurcko": str,
                "MolFormula": str,
                "AtomBondCounts": str,
                "DegreeVector": str,
                "Mesomer": str,
                "HetAtomTautomer": str,
                "HetAtomProtomer": str,
                "RedoxPair": str,
                "Regioisomer": str,
                "NetCharge": str,
                "SmallWorldIndexBR": str,
                "SmallWorldIndexBRL": str,
                "ArthorSubstructureOrder": str,
                "HetAtomTautomerv2": str,
                "inchi": str,
                "CanonicalSmiles": str,
            },
        )


class ContainsRel(Relationship):
    def __init__(self, *args, **kwargs):
        super().__init__(
            neo4j_label="CONTAINS",
            description="Connects a chemical compound to a chemical fragment",
            start=ChemicalCompound,
            end=ChemicalFragment,
        )