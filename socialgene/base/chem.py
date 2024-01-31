from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolHash
from rdkit import Chem
from rdkit.Chem import rdBase
from socialgene.utils.logging import log
from rdkit.Chem import Descriptors


class ChemicalCompound:
    def __init__(self, compound):
        self.mol = None
        self.parse_compound(compound)
        # self.uid = AllChem.GetMorganFingerprintAsBitVect(self.mol, 2)

    def parse_compound(self, input):
        if isinstance(input, Chem.rdchem.Mol):
            self.mol = input
        elif isinstance(input, str):
            method_list = [Chem.MolFromSmiles, Chem.MolFromInchi, Chem.MolFromMolFile]
            for method in method_list:
                log.info(f"Trying to parse compound with {method.__name__}")
                try:
                    temp = method(input)
                    if isinstance(temp, Chem.rdchem.Mol):
                        self.mol = temp
                        log.info(f"Successfully parsed compound with {method.__name__}")
                        break
                except:
                    continue
        if not isinstance(self.mol, Chem.rdchem.Mol):
            raise ValueError("Wasn't able to parse the compound")

    @property
    def hash_dict(self):
        return {
            k: rdMolHash.MolHash(self.mol, v)
            for k, v in rdMolHash.HashFunction.names.items()
        }

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
            except:
                out_dict[i] = None
        return out_dict

    @property
    def fragments(self):
        out_dict = {}
        for i in (i for i, j in Descriptors._descList if i.startswith("fr_")):
            try:
                out_dict[i] = getattr(Descriptors, i)(self.mol)
            except:
                out_dict[i] = None
        return out_dict
