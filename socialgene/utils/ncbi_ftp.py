from ftplib import FTP
from io import BytesIO

from socialgene.utils.logging import log


class NcbiFtp:
    ftp = FTP("ftp.ncbi.nlm.nih.gov", "anonymous", "cmclark8@wisc.edu", timeout=60)

    def __init__(self):
        log.debug(
            f'Connected to "{self.ftp.host}" on {self.ftp.port};\n\n NCBI FTP info: {self.ftp.welcome.replace("220", "")}'
        )
        self.ftp.set_pasv(True)
        self.assembly_paths = []

    def _connect_to_ftp(self):
        log.debug("Connecting to NCBI FTP")
        self.ftp.login()

    def get_assembly_dirpaths(self):
        with BytesIO() as r:
            log.debug("Beginning download of assembly_summary_refseq.txt")
            self.ftp.retrbinary(
                "RETR /genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt", r.write
            )
            self.assembly_paths = []
            for i in r.getvalue().decode("utf-8").split("\n"):
                if i.startswith("#"):
                    pass
                elif i == "":
                    pass
                else:
                    self.assembly_paths.append(i.split("\t")[19])
        self.assembly_paths = [
            i.replace("ftp://ftp.ncbi.nlm.nih.gov", "") for i in self.assembly_paths
        ]
        self.assembly_paths = [
            i.replace("https://ftp.ncbi.nlm.nih.gov", "") for i in self.assembly_paths
        ]
        self.assembly_paths = [i.strip("/") for i in self.assembly_paths]
        log.debug("Finished download of assembly_summary_refseq.txt")
        log.debug(f"Found {len(self.assembly_paths)} assembly urls")
