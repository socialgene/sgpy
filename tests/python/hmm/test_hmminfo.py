import os
from pathlib import Path
import tempfile
from socialgene.hmm.hmminfo import HmmInfo
from socialgene.utils.logging import log

# FIXTURE_DIR = os.path.dirname(os.path.realpath(__file__))
# FIXTURE_DIR = os.path.dirname(FIXTURE_DIR)
# FIXTURE_DIR = os.path.join(FIXTURE_DIR, "data")
FIXTURE_DIR = "/home/chase/Documents/github/kwan_lab/socialgene/sgpy/tests/python/data"
hmm_metadata = os.path.join(FIXTURE_DIR, "hmm_metadata.tsv.gz")


def test_read():
    expected = {
        (
            "amrfinder",
            "amrfinder/156720500_bleO-NCBIFAM.HMM_so",
            "156720500_bleO-NCBIFAM",
            "NF000027.1",
            "NCBIFAM: bleomycin binding protein",
            "Tue Aug 27 12:35:11 2019",
            "gTBtzgszynDnv7UBSOAna-E_VpphlMxy",
            "134",
            "",
            "",
        ),
        (
            "antismash",
            "antismash/antismash/detection/genefunctions/data/smcogs.hmm_so",
            "SMCOG1000:ABC_transporter_ATP-binding_protein",
            "",
            "",
            "Mon Aug 16 11:49:58 2010",
            "rCUJzpN0AqWVMT-59zp3e3IpFQnYj2DB",
            "220",
            "genefunctions",
            "",
        ),
        (
            "resfams",
            "resfams/resfams.hmm_so",
            "AAC3-I",
            "RF0003",
            "Aminoglycoside Acetyltransferase (AAC3-I) [ARO:3001205]",
            "Mon Jun 24 15:47:02 2013",
            "MPBuSHykeoJK_ikb1-GjJuXQB30UcTd7",
            "154",
            "",
            "",
        ),
        (
            "pfam",
            "pfam/Pfam-A.hmm_so",
            "1-cysPrx_C",
            "PF10417.12",
            "C-terminal domain of 1-Cys peroxiredoxin",
            "Thu Nov  4 19:20:02 2021",
            "cEWBLaPObSKpbnmbPlxstrBDrkoCy7J6",
            "40",
            "",
            "",
        ),
        (
            "amrfinder",
            "amrfinder/154989_fosA-NCBIFAM.HMM_so",
            "154989_fosA-NCBIFAM",
            "NF000026.1",
            "NCBIFAM: FosA/FosA2 family fosfomycin resistance glutathione transferase",
            "Tue Aug 27 12:35:11 2019",
            "KJuPYmP-uSiSU-yKk1AqNRBBROGFXO4y",
            "141",
            "",
            "",
        ),
        (
            "antismash",
            "antismash/antismash/detection/genefunctions/data/smcogs.hmm_so",
            "SMCOG1001:short-chain_dehydrogenase/reductase_SDR",
            "",
            "",
            "Mon Aug 16 11:50:11 2010",
            "qtnf2MTnCWbo8bkRPzQ2ifkVfyaaEi6B",
            "246",
            "genefunctions",
            "",
        ),
        (
            "tigrfam",
            "tigrfam/TIGRFAMs_15.0_HMM.hmm_so",
            "TIGR00002",
            "TIGR00002",
            "S16: ribosomal protein bS16",
            "Sun Apr 11 21:57:51 2010",
            "Oc7q-XLdzv60_fuySzvRbS52BwNwJQ31",
            "78",
            "",
            "",
        ),
        (
            "tigrfam",
            "tigrfam/TIGRFAMs_15.0_HMM.hmm_so",
            "TIGR00003",
            "TIGR00003",
            "TIGR00003: copper ion binding protein",
            "Sun Apr 11 21:57:51 2010",
            "f_R2niSD_YOB1igPvoJZxZdDMa9_olgq",
            "66",
            "",
            "",
        ),
        (
            "pfam",
            "pfam/Pfam-A.hmm_so",
            "120_Rick_ant",
            "PF12574.11",
            "120 KDa Rickettsia surface antigen",
            "Tue Oct 12 02:07:11 2021",
            "YaH4-hg-bGWD7YANUul1HZbSlVdCRptX",
            "238",
            "",
            "",
        ),
        (
            "resfams",
            "resfams/resfams.hmm_so",
            "16S_rRNA_methyltrans",
            "RF0001",
            "16S ribosomal RNA methyltransferase [ARO:3000857]",
            "Mon Jun 24 15:56:03 2013",
            "SHfkTFMKLE4JFTkGbd8PWf5-dxvcxSbE",
            "251",
            "",
            "",
        ),
        (
            "amrfinder",
            "amrfinder/1567214_ble-NCBIFAM.HMM_so",
            "1567214_ble-NCBIFAM",
            "NF000028.1",
            "NCBIFAM: BLMA family bleomycin binding protein",
            "Tue Aug 27 12:35:11 2019",
            "Wq6vIHD8UjHbXXtCT5Eh2ELziqPYmx87",
            "122",
            "",
            "",
        ),
        (
            "pfam",
            "pfam/Pfam-A.hmm_so",
            "12TM_1",
            "PF09847.12",
            "Membrane protein of 12 TMs",
            "Wed Oct 20 15:13:22 2021",
            "mXr20sOTu_bX57pxUIiMk7Pa4xfsUme3",
            "449",
            "",
            "",
        ),
        (
            "resfams",
            "resfams/resfams.hmm_so",
            "AAC3",
            "RF0002",
            "Aminoglycoside Acetyltransferase (AAC3) [ARO:3000322]",
            "Mon Jun 24 15:56:52 2013",
            "d0iuEX_bl-bMxcd9rLO1_D2akDhTgdrX",
            "292",
            "",
            "",
        ),
        (
            "tigrfam",
            "tigrfam/TIGRFAMs_15.0_HMM.hmm_so",
            "TIGR00001",
            "TIGR00001",
            "rpmI_bact: ribosomal protein bL35",
            "Sun Apr 11 21:57:51 2010",
            "MyTtcT35uXBcA2k-osQXnhmjCf7qzLUJ",
            "63",
            "",
            "",
        ),
        (
            "antismash",
            "antismash/antismash/detection/genefunctions/data/smcogs.hmm_so",
            "SMCOG1002:AMP-dependent_synthetase_and_ligase",
            "",
            "",
            "Mon Aug 16 11:51:26 2010",
            "Fpg1dfiXGoZWiF5Hg80pRUvsabYvVlmJ",
            "498",
            "genefunctions",
            "",
        ),
    }
    a = HmmInfo(hmm_metadata)
    a.read_tsv()
    assert expected == a.all_hmms_data
    assert 15 == len(a.all_hmms_data)


def test_write_nr_hmm_nodes():
    expected = [
        "Fpg1dfiXGoZWiF5Hg80pRUvsabYvVlmJ\t498\n",
        "KJuPYmP-uSiSU-yKk1AqNRBBROGFXO4y\t141\n",
        "MPBuSHykeoJK_ikb1-GjJuXQB30UcTd7\t154\n",
        "MyTtcT35uXBcA2k-osQXnhmjCf7qzLUJ\t63\n",
        "Oc7q-XLdzv60_fuySzvRbS52BwNwJQ31\t78\n",
        "SHfkTFMKLE4JFTkGbd8PWf5-dxvcxSbE\t251\n",
        "Wq6vIHD8UjHbXXtCT5Eh2ELziqPYmx87\t122\n",
        "YaH4-hg-bGWD7YANUul1HZbSlVdCRptX\t238\n",
        "cEWBLaPObSKpbnmbPlxstrBDrkoCy7J6\t40\n",
        "d0iuEX_bl-bMxcd9rLO1_D2akDhTgdrX\t292\n",
        "f_R2niSD_YOB1igPvoJZxZdDMa9_olgq\t66\n",
        "gTBtzgszynDnv7UBSOAna-E_VpphlMxy\t134\n",
        "mXr20sOTu_bX57pxUIiMk7Pa4xfsUme3\t449\n",
        "qtnf2MTnCWbo8bkRPzQ2ifkVfyaaEi6B\t246\n",
        "rCUJzpN0AqWVMT-59zp3e3IpFQnYj2DB\t220\n",
    ]
    a = HmmInfo(hmm_metadata)
    a.read_tsv()
    with tempfile.NamedTemporaryFile() as fp:
        a.write_nr_hmm_nodes(fp.name)
        with open(fp.name, "r") as h:
            z = h.readlines()
    z.sort()
    assert expected == z


def test_write_hmm_source_nodes():
    expected_fnames = [
        "amrfinder_hmm_source",
        "antismash_hmm_source",
        "pfam_hmm_source",
        "resfams_hmm_source",
        "tigrfam_hmm_source",
    ]
    expected_file_contents = [
        "amrfinder/154989_fosA-NCBIFAM.HMM_so\t154989_fosA-NCBIFAM\tNF000026.1\tNCBIFAM: FosA/FosA2 family fosfomycin resistance glutathione transferase\tTue Aug 27 12:35:11 2019\tKJuPYmP-uSiSU-yKk1AqNRBBROGFXO4y\t141\t\t\n",
        "amrfinder/156720500_bleO-NCBIFAM.HMM_so\t156720500_bleO-NCBIFAM\tNF000027.1\tNCBIFAM: bleomycin binding protein\tTue Aug 27 12:35:11 2019\tgTBtzgszynDnv7UBSOAna-E_VpphlMxy\t134\t\t\n",
        "amrfinder/1567214_ble-NCBIFAM.HMM_so\t1567214_ble-NCBIFAM\tNF000028.1\tNCBIFAM: BLMA family bleomycin binding protein\tTue Aug 27 12:35:11 2019\tWq6vIHD8UjHbXXtCT5Eh2ELziqPYmx87\t122\t\t\n",
        "antismash/antismash/detection/genefunctions/data/smcogs.hmm_so\tSMCOG1000:ABC_transporter_ATP-binding_protein\t\t\tMon Aug 16 11:49:58 2010\trCUJzpN0AqWVMT-59zp3e3IpFQnYj2DB\t220\tgenefunctions\t\n",
        "antismash/antismash/detection/genefunctions/data/smcogs.hmm_so\tSMCOG1001:short-chain_dehydrogenase/reductase_SDR\t\t\tMon Aug 16 11:50:11 2010\tqtnf2MTnCWbo8bkRPzQ2ifkVfyaaEi6B\t246\tgenefunctions\t\n",
        "antismash/antismash/detection/genefunctions/data/smcogs.hmm_so\tSMCOG1002:AMP-dependent_synthetase_and_ligase\t\t\tMon Aug 16 11:51:26 2010\tFpg1dfiXGoZWiF5Hg80pRUvsabYvVlmJ\t498\tgenefunctions\t\n",
        "pfam/Pfam-A.hmm_so\t1-cysPrx_C\tPF10417.12\tC-terminal domain of 1-Cys peroxiredoxin\tThu Nov  4 19:20:02 2021\tcEWBLaPObSKpbnmbPlxstrBDrkoCy7J6\t40\t\t\n",
        "pfam/Pfam-A.hmm_so\t120_Rick_ant\tPF12574.11\t120 KDa Rickettsia surface antigen\tTue Oct 12 02:07:11 2021\tYaH4-hg-bGWD7YANUul1HZbSlVdCRptX\t238\t\t\n",
        "pfam/Pfam-A.hmm_so\t12TM_1\tPF09847.12\tMembrane protein of 12 TMs\tWed Oct 20 15:13:22 2021\tmXr20sOTu_bX57pxUIiMk7Pa4xfsUme3\t449\t\t\n",
        "resfams/resfams.hmm_so\t16S_rRNA_methyltrans\tRF0001\t16S ribosomal RNA methyltransferase [ARO:3000857]\tMon Jun 24 15:56:03 2013\tSHfkTFMKLE4JFTkGbd8PWf5-dxvcxSbE\t251\t\t\n",
        "resfams/resfams.hmm_so\tAAC3\tRF0002\tAminoglycoside Acetyltransferase (AAC3) [ARO:3000322]\tMon Jun 24 15:56:52 2013\td0iuEX_bl-bMxcd9rLO1_D2akDhTgdrX\t292\t\t\n",
        "resfams/resfams.hmm_so\tAAC3-I\tRF0003\tAminoglycoside Acetyltransferase (AAC3-I) [ARO:3001205]\tMon Jun 24 15:47:02 2013\tMPBuSHykeoJK_ikb1-GjJuXQB30UcTd7\t154\t\t\n",
        "tigrfam/TIGRFAMs_15.0_HMM.hmm_so\tTIGR00001\tTIGR00001\trpmI_bact: ribosomal protein bL35\tSun Apr 11 21:57:51 2010\tMyTtcT35uXBcA2k-osQXnhmjCf7qzLUJ\t63\t\t\n",
        "tigrfam/TIGRFAMs_15.0_HMM.hmm_so\tTIGR00002\tTIGR00002\tS16: ribosomal protein bS16\tSun Apr 11 21:57:51 2010\tOc7q-XLdzv60_fuySzvRbS52BwNwJQ31\t78\t\t\n",
        "tigrfam/TIGRFAMs_15.0_HMM.hmm_so\tTIGR00003\tTIGR00003\tTIGR00003: copper ion binding protein\tSun Apr 11 21:57:51 2010\tf_R2niSD_YOB1igPvoJZxZdDMa9_olgq\t66\t\t\n",
    ]
    a = HmmInfo(hmm_metadata)
    a.read_tsv()
    with tempfile.TemporaryDirectory() as fp:
        a.write_hmm_source_nodes(fp)
        pz = Path(fp)
        contents = list(pz.glob("*"))
        fnames = [i.stem for i in contents]
        fnames.sort()
        # assert expected_fnames == fnames
        holdit = []
        for fname in fnames:
            with open({i.stem: i for i in contents}[fname], "r") as h:
                holdit.extend(h.readlines())
        holdit.sort()
        assert expected_file_contents == holdit


def bro():
    a = HmmInfo(hmm_metadata)
    a.read_tsv()
    with tempfile.TemporaryDirectory() as fp:
        a.write_hmm_source_nodes(fp)
        pz = Path(fp)
        contents = list(pz.glob("*"))
        fnames = [i.stem for i in contents]
        fnames.sort()
        # assert expected_fnames == fnames
        holdit = []
        for fname in fnames:
            with open({i.stem: i for i in contents}[fname], "r") as h:
                holdit.extend(h.readlines())
        return holdit
