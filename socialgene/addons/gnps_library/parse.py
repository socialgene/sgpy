from socialgene.addons.gnps_library.nr import GnpsLibrarySpectrumNode


class GNPSLibrarySpectrum:
    def __init__(self, specnet: dict):
        self.properties = {}

    def from_specnet(self, specnet: dict):
        self.properties["uid"] = str(specnet["spectrumid"])
        self.properties["compound_name"] = str(specnet["compound_name"])
        self.properties["compound_source"] = str(specnet["compound_source"])
        self.properties["pi"] = str(specnet["pi"])
        self.properties["data_collector"] = str(specnet["data_collector"])
        self.properties["adduct"] = str(specnet["adduct"])
        self.properties["precursor_mz"] = float(specnet["precursor_mz"])
        self.properties["exactmass"] = float(specnet["exactmass"])
        self.properties["charge"] = int(specnet["charge"])
        self.properties["cas_number"] = str(specnet["cas_number"])
        self.properties["pubmed_id"] = str(specnet["pubmed_id"])
        self.properties["smiles"] = str(specnet["smiles"])
        self.properties["inchi"] = str(specnet["inchi"])
        self.properties["inchi_aux"] = str(specnet["inchi_aux"])
        self.properties["library_class"] = str(specnet["library_class"])
        self.properties["ionmode"] = str(specnet["ionmode"])
        self.properties["libraryqualitystring"] = str(specnet["libraryqualitystring"])
        self.properties["mqscore"] = float(specnet["mqscore"])
        self.properties["tic_query"] = float(specnet["tic_query"])
        self.properties["rt_query"] = float(specnet["rt_query"])
        self.properties["mzerrorppm"] = float(specnet["mzerrorppm"])
        self.properties["sharedpeaks"] = int(specnet["sharedpeaks"])
        self.properties["massdiff"] = float(specnet["massdiff"])
        self.properties["libmz"] = float(specnet["libmz"])
        self.properties["specmz"] = float(specnet["specmz"])
        self.properties["speccharge"] = int(specnet["speccharge"])
        self.properties["moleculeexplorerdatasets"] = str(
            specnet["moleculeexplorerdatasets"]
        )
        self.properties["moleculeexplorerfiles"] = str(specnet["moleculeexplorerfiles"])
        self.properties["molecular_formula"] = str(specnet["molecular_formula"])
        self.properties["inchikey"] = str(specnet["inchikey"])
        self.properties["inchikey_planar"] = str(specnet["inchikey_planar"])

    @property
    def node(self):
        if self.properties:
            return GnpsLibrarySpectrumNode(properties=self.properties)
