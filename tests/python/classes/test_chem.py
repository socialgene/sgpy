import pytest
from socialgene.base.chem import ChemicalCompound
from rdkit import Chem


EXPECTED_FRAGMENTS = (
    "fr_Al_COO",
    "fr_Al_OH",
    "fr_Al_OH_noTert",
    "fr_ArN",
    "fr_Ar_COO",
    "fr_Ar_N",
    "fr_Ar_NH",
    "fr_Ar_OH",
    "fr_COO",
    "fr_COO2",
    "fr_C_O",
    "fr_C_O_noCOO",
    "fr_C_S",
    "fr_HOCCN",
    "fr_Imine",
    "fr_NH0",
    "fr_NH1",
    "fr_NH2",
    "fr_N_O",
    "fr_Ndealkylation1",
    "fr_Ndealkylation2",
    "fr_Nhpyrrole",
    "fr_SH",
    "fr_aldehyde",
    "fr_alkyl_carbamate",
    "fr_alkyl_halide",
    "fr_allylic_oxid",
    "fr_amide",
    "fr_amidine",
    "fr_aniline",
    "fr_aryl_methyl",
    "fr_azide",
    "fr_azo",
    "fr_barbitur",
    "fr_benzene",
    "fr_benzodiazepine",
    "fr_bicyclic",
    "fr_diazo",
    "fr_dihydropyridine",
    "fr_epoxide",
    "fr_ester",
    "fr_ether",
    "fr_furan",
    "fr_guanido",
    "fr_halogen",
    "fr_hdrzine",
    "fr_hdrzone",
    "fr_imidazole",
    "fr_imide",
    "fr_isocyan",
    "fr_isothiocyan",
    "fr_ketone",
    "fr_ketone_Topliss",
    "fr_lactam",
    "fr_lactone",
    "fr_methoxy",
    "fr_morpholine",
    "fr_nitrile",
    "fr_nitro",
    "fr_nitro_arom",
    "fr_nitro_arom_nonortho",
    "fr_nitroso",
    "fr_oxazole",
    "fr_oxime",
    "fr_para_hydroxylation",
    "fr_phenol",
    "fr_phenol_noOrthoHbond",
    "fr_phos_acid",
    "fr_phos_ester",
    "fr_piperdine",
    "fr_piperzine",
    "fr_priamide",
    "fr_prisulfonamd",
    "fr_pyridine",
    "fr_quatN",
    "fr_sulfide",
    "fr_sulfonamd",
    "fr_sulfone",
    "fr_term_acetylene",
    "fr_tetrazole",
    "fr_thiazole",
    "fr_thiocyan",
    "fr_thiophene",
    "fr_unbrch_alkane",
    "fr_urea",
)


aspirin_hash_dict = {
    "AnonymousGraph": "**(*)**1*****1*(*)*",
    "ElementGraph": "CC(O)OC1CCCCC1C(O)O",
    "CanonicalSmiles": "CC(=O)Oc1ccccc1C(=O)O",
    "MurckoScaffold": "c1ccccc1",
    "ExtendedMurcko": "*c1ccccc1*",
    "MolFormula": "C9H8O4",
    "AtomBondCounts": "13,13",
    "DegreeVector": "0,4,5,4",
    "Mesomer": "C[C]([O])O[C]1[CH][CH][CH][CH][C]1[C]([O])O_0",
    "HetAtomTautomer": "C[C]([O])O[C]1[CH][CH][CH][CH][C]1[C]([O])[O]_1_0",
    "HetAtomProtomer": "C[C]([O])O[C]1[CH][CH][CH][CH][C]1[C]([O])[O]_1",
    "RedoxPair": "C[C]([O])O[C]1[CH][CH][CH][CH][C]1[C]([O])O",
    "Regioisomer": "*C(=O)O.*OC(C)=O.c1ccccc1",
    "NetCharge": "0",
    "SmallWorldIndexBR": "B13R1",
    "SmallWorldIndexBRL": "B13R1L5",
    "ArthorSubstructureOrder": "000d000d0100090004000056000000",
    "HetAtomTautomerv2": "[C]:[C](:[O]):[O]:[C]1:[C]:[C]:[C]:[C]:[C]:1:[C](:[O]):[O]_8_0",
}

aspirin_base_properties = {
    "MolWt": 180.15899999999996,
    "HeavyAtomMolWt": 172.09499999999997,
    "ExactMolWt": 180.042258736,
    "NumValenceElectrons": 68,
    "NumRadicalElectrons": 0,
    "HeavyAtomCount": 13,
    "NumAliphaticCarbocycles": 0,
    "NumAliphaticHeterocycles": 0,
    "NumAliphaticRings": 0,
    "NumAromaticCarbocycles": 1,
    "NumAromaticHeterocycles": 0,
    "NumAromaticRings": 1,
    "NumHAcceptors": 3,
    "NumHDonors": 1,
    "NumHeteroatoms": 4,
    "NumRotatableBonds": 2,
    "NumSaturatedCarbocycles": 0,
    "NumSaturatedHeterocycles": 0,
    "NumSaturatedRings": 0,
    "RingCount": 1,
    "MolLogP": 1.3101,
    "MolMR": 44.71030000000002,
}

aspirin_frags = {
    "fr_Al_COO": 0,
    "fr_Al_OH": 0,
    "fr_Al_OH_noTert": 0,
    "fr_ArN": 0,
    "fr_Ar_COO": 1,
    "fr_Ar_N": 0,
    "fr_Ar_NH": 0,
    "fr_Ar_OH": 0,
    "fr_COO": 1,
    "fr_COO2": 1,
    "fr_C_O": 2,
    "fr_C_O_noCOO": 1,
    "fr_C_S": 0,
    "fr_HOCCN": 0,
    "fr_Imine": 0,
    "fr_NH0": 0,
    "fr_NH1": 0,
    "fr_NH2": 0,
    "fr_N_O": 0,
    "fr_Ndealkylation1": 0,
    "fr_Ndealkylation2": 0,
    "fr_Nhpyrrole": 0,
    "fr_SH": 0,
    "fr_aldehyde": 0,
    "fr_alkyl_carbamate": 0,
    "fr_alkyl_halide": 0,
    "fr_allylic_oxid": 0,
    "fr_amide": 0,
    "fr_amidine": 0,
    "fr_aniline": 0,
    "fr_aryl_methyl": 0,
    "fr_azide": 0,
    "fr_azo": 0,
    "fr_barbitur": 0,
    "fr_benzene": 1,
    "fr_benzodiazepine": 0,
    "fr_bicyclic": 0,
    "fr_diazo": 0,
    "fr_dihydropyridine": 0,
    "fr_epoxide": 0,
    "fr_ester": 1,
    "fr_ether": 1,
    "fr_furan": 0,
    "fr_guanido": 0,
    "fr_halogen": 0,
    "fr_hdrzine": 0,
    "fr_hdrzone": 0,
    "fr_imidazole": 0,
    "fr_imide": 0,
    "fr_isocyan": 0,
    "fr_isothiocyan": 0,
    "fr_ketone": 0,
    "fr_ketone_Topliss": 0,
    "fr_lactam": 0,
    "fr_lactone": 0,
    "fr_methoxy": 0,
    "fr_morpholine": 0,
    "fr_nitrile": 0,
    "fr_nitro": 0,
    "fr_nitro_arom": 0,
    "fr_nitro_arom_nonortho": 0,
    "fr_nitroso": 0,
    "fr_oxazole": 0,
    "fr_oxime": 0,
    "fr_para_hydroxylation": 1,
    "fr_phenol": 0,
    "fr_phenol_noOrthoHbond": 0,
    "fr_phos_acid": 0,
    "fr_phos_ester": 0,
    "fr_piperdine": 0,
    "fr_piperzine": 0,
    "fr_priamide": 0,
    "fr_prisulfonamd": 0,
    "fr_pyridine": 0,
    "fr_quatN": 0,
    "fr_sulfide": 0,
    "fr_sulfonamd": 0,
    "fr_sulfone": 0,
    "fr_term_acetylene": 0,
    "fr_tetrazole": 0,
    "fr_thiazole": 0,
    "fr_thiocyan": 0,
    "fr_thiophene": 0,
    "fr_unbrch_alkane": 0,
    "fr_urea": 0,
}


def test_ChemicalCompound_smiles():
    compound = ChemicalCompound("CC(=O)OC1=CC=CC=C1C(=O)O")
    assert isinstance(compound.mol, Chem.rdchem.Mol)
    assert compound.hash_dict == aspirin_hash_dict
    assert pytest.approx(compound.base_properties) == aspirin_base_properties
    assert compound.fragments == aspirin_frags


def test_ChemicalCompound_rdchemmolinchi():
    compound = ChemicalCompound(
        "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    )
    assert isinstance(compound.mol, Chem.rdchem.Mol)
    assert compound.hash_dict == aspirin_hash_dict
    assert aspirin_base_properties == pytest.approx(compound.base_properties)
    assert compound.fragments == aspirin_frags


def test_ChemicalCompound_rdchemmol():
    inp = Chem.MolFromInchi(
        "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    )
    compound = ChemicalCompound(inp)
    assert isinstance(compound.mol, Chem.rdchem.Mol)
    assert compound.hash_dict == aspirin_hash_dict
    assert aspirin_base_properties == pytest.approx(compound.base_properties)
    assert compound.fragments == aspirin_frags


def test_ChemicalCompound_invalid_format():
    with pytest.raises(ValueError):
        ChemicalCompound("1")


def test_ChemicalCompound_smiles2():
    compound = ChemicalCompound(
        "CCCCCCCCCC(=O)NC(CC1=CNC2=CC=CC=C21)C(=O)NC(CC(=O)N)C(=O)NC(CC(=O)O)C(=O)NC3C(OC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)C(C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"
    )
    assert isinstance(compound.mol, Chem.rdchem.Mol)
    assert compound.hash_dict == {
        "AnonymousGraph": "**********(*)**(**1***2*****12)*(*)**(**(*)*)*(*)**(**(*)*)*(*)**1*(*)***(*)**(****)*(*)**(**(*)*)*(*)**(*)*(*)**(**(*)*)*(*)***(*)**(**)*(*)**(*(*)**(*)*)*(*)**(**(*)*2*****2*)*(*)**1*",
        "ElementGraph": "CCCCCCCCCC(O)NC(CC1CNC2CCCCC12)C(O)NC(CC(N)O)C(O)NC(CC(O)O)C(O)NC1C(O)NCC(O)NC(CCCN)C(O)NC(CC(O)O)C(O)NC(C)C(O)NC(CC(O)O)C(O)NCC(O)NC(CO)C(O)NC(C(C)CC(O)O)C(O)NC(CC(O)C2CCCCC2N)C(O)OC1C",
        "CanonicalSmiles": "CCCCCCCCCC(=O)NC(Cc1c[nH]c2ccccc12)C(=O)NC(CC(N)=O)C(=O)NC(CC(=O)O)C(=O)NC1C(=O)NCC(=O)NC(CCCN)C(=O)NC(CC(=O)O)C(=O)NC(C)C(=O)NC(CC(=O)O)C(=O)NCC(=O)NC(CO)C(=O)NC(C(C)CC(=O)O)C(=O)NC(CC(=O)c2ccccc2N)C(=O)OC1C",
        "MurckoScaffold": "c1ccc(CCC2COCC(NCCNCCNCCCc3c[nH]c4ccccc34)CNCCNCCNCCNCCNCCNCCNCCNCCN2)cc1",
        "ExtendedMurcko": "*c1ccccc1C(=*)CC1NC(=*)C(*)NC(=*)C(*)NC(=*)CNC(=*)C(*)NC(=*)C(*)NC(=*)C(*)NC(=*)C(*)NC(=*)CNC(=*)C(NC(=*)C(*)NC(=*)C(*)NC(=*)C(*)Cc2c[nH]c3ccccc23)C(*)OC1=*",
        "MolFormula": "C72H101N17O26",
        "AtomBondCounts": "115,118",
        "DegreeVector": "0,38,45,32",
        "Mesomer": "CCCCCCCCC[C]([O])NC(C[C]1[CH]N[C]2[CH][CH][CH][CH][C]12)[C]([O])NC(C[C](N)[O])[C]([O])NC(C[C]([O])O)[C]([O])NC1[C]([O])NC[C]([O])NC(CCCN)[C]([O])NC(C[C]([O])O)[C]([O])NC(C)[C]([O])NC(C[C]([O])O)[C]([O])NC[C]([O])NC(CO)[C]([O])NC(C(C)C[C]([O])O)[C]([O])NC(C[C]([O])[C]2[CH][CH][CH][CH][C]2N)[C]([O])OC1C_0",
        "HetAtomTautomer": "CCCCCCCCC[C]([O])[N]C(C[C]1[CH][N][C]2[CH][CH][CH][CH][C]12)[C]([O])[N]C(C[C]([N])[O])[C]([O])[N]C(C[C]([O])[O])[C]([O])[N]C1[C]([O])[N]C[C]([O])[N]C(CCC[N])[C]([O])[N]C(C[C]([O])[O])[C]([O])[N]C(C)[C]([O])[N]C(C[C]([O])[O])[C]([O])[N]C[C]([O])[N]C(C[O])[C]([O])[N]C(C(C)C[C]([O])[O])[C]([O])[N]C(C[C]([O])[C]2[CH][CH][CH][CH][C]2[N])[C]([O])OC1C_25_0",
        "HetAtomProtomer": "CCCCCCCCC[C]([O])[N]C(C[C]1[CH][N][C]2[CH][CH][CH][CH][C]12)[C]([O])[N]C(C[C]([N])[O])[C]([O])[N]C(C[C]([O])[O])[C]([O])[N]C1[C]([O])[N]C[C]([O])[N]C(CCC[N])[C]([O])[N]C(C[C]([O])[O])[C]([O])[N]C(C)[C]([O])[N]C(C[C]([O])[O])[C]([O])[N]C[C]([O])[N]C(C[O])[C]([O])[N]C(C(C)C[C]([O])[O])[C]([O])[N]C(C[C]([O])[C]2[CH][CH][CH][CH][C]2[N])[C]([O])OC1C_25",
        "RedoxPair": "CCCCCCCCC[C]([O])NC(C[C]1[CH]N[C]2[CH][CH][CH][CH][C]12)[C]([O])NC(C[C](N)[O])[C]([O])NC(C[C]([O])O)[C]([O])NC1[C]([O])NC[C]([O])NC(CCCN)[C]([O])NC(C[C]([O])O)[C]([O])NC(C)[C]([O])NC(C[C]([O])O)[C]([O])NC[C]([O])NC(CO)[C]([O])NC(C(C)C[C]([O])O)[C]([O])NC(C[C]([O])[C]2[CH][CH][CH][CH][C]2N)[C]([O])OC1C",
        "Regioisomer": "*C.*C.*C.*C(C)CC(=O)O.*CC(*)=O.*CC(=O)O.*CC(=O)O.*CCC.*CCC(=O)N*.*N.*N.*NC(=O)CCC(=O)O.*NC(=O)CCC(N)=O.*NC(=O)CCCCCCCCC.*O.O=C1CCOC(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)CN1.c1ccc2[nH]ccc2c1.c1ccccc1",
        "NetCharge": "0",
        "SmallWorldIndexBR": "B118R4",
        "SmallWorldIndexBRL": "B118R4L45",
        "ArthorSubstructureOrder": "00730076010048002b0002f7000000",
        "HetAtomTautomerv2": "[CH3]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[C]:[C](:[O]):[N]:[C](-[CH2]-[C]1:[C]:[N]:[C]2:[C]:[C]:[C]:[C]:[C]:1:2):[C](:[O]):[N]:[C](-[C]:[C](:[N]):[O]):[C](:[O]):[N]:[C](-[C]:[C](:[O]):[O]):[C](:[O]):[N]:[C]1:[C](:[O]):[N]:[C]:[C](:[O]):[N]:[C](-[CH2]-[CH2]-[CH2]-[NH2]):[C](:[O]):[N]:[C](-[C]:[C](:[O]):[O]):[C](:[O]):[N]:[C](-[CH3]):[C](:[O]):[N]:[C](-[C]:[C](:[O]):[O]):[C](:[O]):[N]:[C]:[C](:[O]):[N]:[C](-[CH2]-[OH]):[C](:[O]):[N]:[C](-[CH](-[CH3])-[C]:[C](:[O]):[O]):[C](:[O]):[N]:[C](-[C]:[C](:[O]):[C]2:[C]:[C]:[C]:[C]:[C]:2:[N]):[C](:[O]):[O]-[CH]-1-[CH3]_60_0",
    }
    assert compound.base_properties == pytest.approx(
        {
            "MolWt": 1620.693000000001,
            "HeavyAtomMolWt": 1518.884999999999,
            "ExactMolWt": 1619.7103663519993,
            "NumValenceElectrons": 630,
            "NumRadicalElectrons": 0,
            "HeavyAtomCount": 115,
            "NumAliphaticCarbocycles": 0,
            "NumAliphaticHeterocycles": 1,
            "NumAliphaticRings": 1,
            "NumAromaticCarbocycles": 2,
            "NumAromaticHeterocycles": 1,
            "NumAromaticRings": 3,
            "NumHAcceptors": 24,
            "NumHDonors": 22,
            "NumHeteroatoms": 43,
            "NumRotatableBonds": 35,
            "NumSaturatedCarbocycles": 0,
            "NumSaturatedHeterocycles": 1,
            "NumSaturatedRings": 1,
            "RingCount": 4,
            "MolLogP": -5.621799999999963,
            "MolMR": 400.26050000000095,
        }
    )
    assert compound.fragments == {
        "fr_Al_COO": 4,
        "fr_Al_OH": 1,
        "fr_Al_OH_noTert": 1,
        "fr_ArN": 1,
        "fr_Ar_COO": 0,
        "fr_Ar_N": 1,
        "fr_Ar_NH": 1,
        "fr_Ar_OH": 0,
        "fr_COO": 4,
        "fr_COO2": 4,
        "fr_C_O": 20,
        "fr_C_O_noCOO": 16,
        "fr_C_S": 0,
        "fr_HOCCN": 0,
        "fr_Imine": 0,
        "fr_NH0": 0,
        "fr_NH1": 14,
        "fr_NH2": 3,
        "fr_N_O": 0,
        "fr_Ndealkylation1": 0,
        "fr_Ndealkylation2": 0,
        "fr_Nhpyrrole": 1,
        "fr_SH": 0,
        "fr_aldehyde": 0,
        "fr_alkyl_carbamate": 0,
        "fr_alkyl_halide": 0,
        "fr_allylic_oxid": 0,
        "fr_amide": 14,
        "fr_amidine": 0,
        "fr_aniline": 1,
        "fr_aryl_methyl": 0,
        "fr_azide": 0,
        "fr_azo": 0,
        "fr_barbitur": 0,
        "fr_benzene": 2,
        "fr_benzodiazepine": 0,
        "fr_bicyclic": 1,
        "fr_diazo": 0,
        "fr_dihydropyridine": 0,
        "fr_epoxide": 0,
        "fr_ester": 1,
        "fr_ether": 1,
        "fr_furan": 0,
        "fr_guanido": 0,
        "fr_halogen": 0,
        "fr_hdrzine": 0,
        "fr_hdrzone": 0,
        "fr_imidazole": 0,
        "fr_imide": 0,
        "fr_isocyan": 0,
        "fr_isothiocyan": 0,
        "fr_ketone": 1,
        "fr_ketone_Topliss": 1,
        "fr_lactam": 0,
        "fr_lactone": 1,
        "fr_methoxy": 0,
        "fr_morpholine": 0,
        "fr_nitrile": 0,
        "fr_nitro": 0,
        "fr_nitro_arom": 0,
        "fr_nitro_arom_nonortho": 0,
        "fr_nitroso": 0,
        "fr_oxazole": 0,
        "fr_oxime": 0,
        "fr_para_hydroxylation": 2,
        "fr_phenol": 0,
        "fr_phenol_noOrthoHbond": 0,
        "fr_phos_acid": 0,
        "fr_phos_ester": 0,
        "fr_piperdine": 0,
        "fr_piperzine": 0,
        "fr_priamide": 1,
        "fr_prisulfonamd": 0,
        "fr_pyridine": 0,
        "fr_quatN": 0,
        "fr_sulfide": 0,
        "fr_sulfonamd": 0,
        "fr_sulfone": 0,
        "fr_term_acetylene": 0,
        "fr_tetrazole": 0,
        "fr_thiazole": 0,
        "fr_thiocyan": 0,
        "fr_thiophene": 0,
        "fr_unbrch_alkane": 5,
        "fr_urea": 0,
    }

