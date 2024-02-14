import pytest
from rdkit import Chem

from socialgene.base.chem import ChemicalCompound

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
    "Ar_COO": 1,
    "COO": 1,
    "COO2": 1,
    "C_O": 2,
    "C_O_noCOO": 1,
    "benzene": 1,
    "ester": 1,
    "ether": 1,
    "para_hydroxylation": 1,
}


def test_ChemicalCompound_smiles():
    compound = ChemicalCompound("CC(=O)OC1=CC=CC=C1C(=O)O")
    assert isinstance(compound.mol, Chem.rdchem.Mol)
    assert compound.hash_dict == aspirin_hash_dict | {
        "inchi": Chem.MolToInchi(compound.mol)
    }
    assert pytest.approx(compound.base_properties) == aspirin_base_properties
    assert compound.fragments == aspirin_frags


def test_ChemicalCompound_rdchemmolinchi():
    compound = ChemicalCompound(
        "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    )
    assert isinstance(compound.mol, Chem.rdchem.Mol)
    assert compound.hash_dict == aspirin_hash_dict | {
        "inchi": Chem.MolToInchi(compound.mol)
    }
    assert aspirin_base_properties == pytest.approx(compound.base_properties)
    assert compound.fragments == aspirin_frags


def test_ChemicalCompound_rdchemmol():
    inp = Chem.MolFromInchi(
        "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    )
    compound = ChemicalCompound(inp)
    assert isinstance(compound.mol, Chem.rdchem.Mol)
    assert compound.hash_dict == aspirin_hash_dict | {
        "inchi": Chem.MolToInchi(compound.mol)
    }
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
        "inchi": "InChI=1S/C72H101N17O26/c1-5-6-7-8-9-10-11-22-53(93)81-44(25-38-31-76-42-20-15-13-17-39(38)42)66(108)84-45(27-52(75)92)67(109)86-48(30-59(102)103)68(110)89-61-37(4)115-72(114)49(26-51(91)40-18-12-14-19-41(40)74)87-71(113)60(35(2)24-56(96)97)88-69(111)50(34-90)82-55(95)32-77-63(105)46(28-57(98)99)83-62(104)36(3)79-65(107)47(29-58(100)101)85-64(106)43(21-16-23-73)80-54(94)33-78-70(61)112/h12-15,17-20,31,35-37,43-50,60-61,76,90H,5-11,16,21-30,32-34,73-74H2,1-4H3,(H2,75,92)(H,77,105)(H,78,112)(H,79,107)(H,80,94)(H,81,93)(H,82,95)(H,83,104)(H,84,108)(H,85,106)(H,86,109)(H,87,113)(H,88,111)(H,89,110)(H,96,97)(H,98,99)(H,100,101)(H,102,103)",
    } | {"inchi": Chem.MolToInchi(compound.mol)}
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
        "Al_COO": 4,
        "Al_OH": 1,
        "Al_OH_noTert": 1,
        "ArN": 1,
        "Ar_N": 1,
        "Ar_NH": 1,
        "COO": 4,
        "COO2": 4,
        "C_O": 20,
        "C_O_noCOO": 16,
        "NH1": 14,
        "NH2": 3,
        "Nhpyrrole": 1,
        "amide": 14,
        "aniline": 1,
        "benzene": 2,
        "bicyclic": 1,
        "ester": 1,
        "ether": 1,
        "ketone": 1,
        "ketone_Topliss": 1,
        "lactone": 1,
        "para_hydroxylation": 2,
        "priamide": 1,
        "unbrch_alkane": 5,
    }
