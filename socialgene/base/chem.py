from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolHash

from socialgene.addons.chemistry.nr import ChemicalCompoundNode
from socialgene.utils.logging import log

# Turn off C++ warnings when in info mode
if log.level != 20:
    RDLogger.DisableLog("rdApp.*")


class ChemicalFragments:
    def __init__(self, **kwargs):
        for i, j in Descriptors._descList:
            if i.startswith("fr_"):
                setattr(self, i, None)

    def add_mol(self, mol):
        for i, j in Descriptors._descList:
            if i.startswith("fr_"):
                setattr(self, i, j(mol))

    def to_dict(self):
        temp = {k: v for k, v in self.__dict__.items() if k.startswith("fr_") and v > 0}
        temp = {i.removeprefix("fr_"): j for i, j in temp.items()}
        return temp


def morgan_fingerprint(rdkitmol, radius=2, nBits=2048):
    return AllChem.GetMorganFingerprintAsBitVect(
        rdkitmol, useChirality=True, radius=radius, nBits=nBits, bitInfo={}
    )


class ChemicalCompound:
    def __init__(self, compound, **kwargs):
        self.mol = None
        self.parse_compound(compound)
        self.morgan = morgan_fingerprint(self.mol, **kwargs)
        self.inchi = Chem.MolToInchi(self.mol)

    def parse_compound(self, input):
        if isinstance(input, Chem.rdchem.Mol):
            self.mol = input
        elif isinstance(input, str):
            method_list = [Chem.MolFromInchi, Chem.MolFromSmiles, Chem.MolFromMolFile]
            for method in method_list:
                log.debug(f"Trying to parse compound with {method.__name__}")
                try:
                    temp = method(input, sanitize=True)
                    if isinstance(temp, Chem.rdchem.Mol):
                        self.mol = temp
                        log.debug(
                            f"Successfully parsed compound with {method.__name__}"
                        )
                        break
                except Exception:
                    continue
        if not isinstance(self.mol, Chem.rdchem.Mol):
            raise ValueError("Wasn't able to parse the compound")

    @property
    def hash_dict(self):
        temp = {
            k: rdMolHash.MolHash(self.mol, v)
            for k, v in rdMolHash.HashFunction.names.items()
        }
        temp["inchi"] = self.inchi
        return temp

    @property
    def base_properties(self):
        out_dict = {}
        for i in (
            "MolWt",
            "HeavyAtomMolWt",
            "ExactMolWt",
            "NumValenceElectrons",
            "NumRadicalElectrons",
            "HeavyAtomCount",
            "NumAliphaticCarbocycles",
            "NumAliphaticHeterocycles",
            "NumAliphaticRings",
            "NumAromaticCarbocycles",
            "NumAromaticHeterocycles",
            "NumAromaticRings",
            "NumHAcceptors",
            "NumHDonors",
            "NumHeteroatoms",
            "NumRotatableBonds",
            "NumSaturatedCarbocycles",
            "NumSaturatedHeterocycles",
            "NumSaturatedRings",
            "RingCount",
            "MolLogP",
            "MolMR",
        ):
            try:
                out_dict[i] = getattr(Descriptors, i)(self.mol)
            except Exception:
                out_dict[i] = None
        return out_dict

    @property
    def fragments(self):
        fr = ChemicalFragments()
        fr.add_mol(self.mol)
        return fr.to_dict()

    @property
    def node(self):
        temp = ChemicalCompoundNode()
        temp.fill_from_dict(self.base_properties | self.hash_dict)
        return temp


class ChemicalCollection:
    def __init__(self) -> None:
        pass

    def compare_morgan(self, a, b):
        return DataStructs.TanimotoSimilarity(self.morgan, self.morgan)
