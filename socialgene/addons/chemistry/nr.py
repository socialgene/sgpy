from socialgene.neo4j.neo4j_element import Node, Relationship


class ChemicalCompoundNode(Node):
    neo4j_label = ["chemical_compound"]
    description = "Represents a chemical compound"
    required_properties = ["inchi", "CanonicalSmiles"]
    property_specification = {
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
    }
    constraints_unique = ["inchi"]


class TanimotoSimilarity(Relationship):
    neo4j_label = "TANIMOTO_SIMILARITY"
    description = "Connects two chemical compounds that are similar"
    start_class = ChemicalCompoundNode
    end_class = ChemicalCompoundNode
    property_specification = {
        "similarity": int,
    }


class McsSimilarity(Relationship):
    neo4j_label = "MCS_SIMILARITY"
    description = "Connects two chemical compounds that have similar maximium common substructures"
    start_class = ChemicalCompoundNode
    end_class = ChemicalCompoundNode
    property_specification = {
        "similarity": int,
    }


class ChemicalSubstructure(ChemicalCompoundNode):
    neo4j_label = ["substructure"]
    description = "Represents a chemical substructure"
    required_properties = ["inchi", "CanonicalSmiles"]
    property_specification = {
        "uid": str,
        "inchi": str,
        "CanonicalSmiles": str,
    }
    constraints_unique = ["inchi", "CanonicalSmiles"]


class ContainsSubstructure(Relationship):
    neo4j_label = "SUBSTRUCTURE"
    description = "Connects a chemical compound to a chemical substructure"
    start_class = ChemicalCompoundNode
    end_class = ChemicalSubstructure
