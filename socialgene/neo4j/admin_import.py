# python dependencies
import os
import pwd
from pathlib import Path

# external dependencies

# internal dependencies
from socialgene.config import env_vars
from socialgene.utils.run_subprocess import run_subprocess
from socialgene.utils.logging import log
from socialgene.neo4j.sg_modules import SocialgeneModules, Neo4jImportData


class Neo4jAdminImport:
    """Socialgene's Neo4jAdminImport class handles building the Docker/Neo4j admin import command line arguments"""

    def __init__(
        self,
        neo4j_top_dir: str,
        cpus: int = 1,
        additional_args: str = "",
        uid: int = None,
        gid: int = None,
        input_sg_modules: list = None,
        hmmlist: list = None,
        dbms_connector_http_listen_address: int = 7474,
        dbms_connector_bolt_listen_address: int = 7687,
        neo4j_version=env_vars["NEO4J_VERSION"],
        *args,
        **kwargs,
    ):
        """Class for building and executing Neo4j admin import

        Args:
            neo4j_top_dir (str): directory containing the necessary dirs and files for Neo4j admin import
            sg_modules (list): see additional options in `socialgene.neo4j.sg_modules` (default is ["base", "hmms", "ncbi_taxonomy"])
            cpus (int, optional): [description]. Defaults to 1.
            additional_args (str, optional): [description]. Defaults to "".
            uid (int, optional): [description]. Defaults to None.
            gid (int, optional): [description]. Defaults to None.
        """
        super().__init__(*args, **kwargs)
        self.input_sg_modules = input_sg_modules
        self.hmmlist = hmmlist
        self.neo4j_top_dir = neo4j_top_dir
        self.cpus = cpus
        self.additional_args = additional_args
        self.node_relationship_argument_list = []
        self.dbms_connector_http_listen_address = dbms_connector_http_listen_address
        self.dbms_connector_bolt_listen_address = dbms_connector_bolt_listen_address
        self.neo4j_version = neo4j_version
        self.filepaths_to_check = []
        # UID/GID for docker call
        if uid is None:
            self.uid = pwd.getpwuid(os.getuid()).pw_uid
        else:
            self.uid = uid
        if gid is None:
            self.gid = pwd.getpwuid(os.getuid()).pw_gid
        else:
            self.gid = gid

    @staticmethod
    def create_neo4j_directories(neo4j_top_dir):
        """Check Neo4j import directories, if they don't exist, create them or Neo4j's Docker container will create them as root user

        Args:
            neo4j_top_dir (str): directory to make subdirectories in, directory should already contain the Neo4j 'import' directory containing data to be imported into the database

        Raises:
            IOError: _description_
            IOError: _description_
        """
        dirs_to_check = ["data", "logs", "plugins"]
        for single_dir in dirs_to_check:
            dir_to_check = os.path.join(neo4j_top_dir, single_dir)
            if not os.path.exists(dir_to_check):
                # create directory if it doesn't exist, or neo4j will get mad
                os.mkdir(dir_to_check)
                log.info(f"Created directory for Neo4j admin import: {dir_to_check}")
        # plugins doesn't have to be empty (and is actually filled in the nextflow workflow)
        dirs_to_check = ["data", "logs"]
        for single_dir in dirs_to_check:
            dir_to_check = os.path.join(neo4j_top_dir, single_dir)
            if len(os.listdir(dir_to_check)) > 0:
                # make sure the directory is empty
                log.exception(
                    f"Directory ({dir_to_check}) must be empty for Neo4j admin import to work"
                )
                raise IOError
        files_to_check = ["import.report"]
        for single_file in files_to_check:
            file_to_check = os.path.join(neo4j_top_dir, single_file)
            if not os.path.exists(file_to_check):
                # create file if it doesn't exist, or neo4j will get mad
                with open(file_to_check, "w"):
                    pass
                log.info(f"Created file for Neo4j admin import: {file_to_check}")
            if os.path.getsize(file_to_check) > 0:
                # make sure the file is empty
                raise IOError(
                    f"File ({file_to_check}) must be empty for Neo4j admin import to work"
                )

    # while building arguments, make a list of data filepaths to check for
    @staticmethod
    def _single_arg_string_builder(
        label,
        header_filename,
        target_subdirectory,
        target_extension,
        type,
    ):
        """Builds the command line arguments to import data into a new Neo4j database

        Args:
            label (str): neo4j admin import "label"
            header_filename (str): neo4j admin import header filename
            target_subdirectory ([type]): subdirectory containing file to import
            extension (str): extension of the data fle for neo4j admin import
            type (str): "relationships" or "nodes"

        Returns:
            list: [first_part_of_arg_string, header_path_string, data_glob_string] want mutable because will check in later step for gz and append if needed
        """
        # chr(92) is a workaround to insert '\\'
        # return f"--{type}={label}=import/neo4j_headers/{header_filename},import/{target_subdirectory}/^.*{chr(92)}.{target_extension}"
        return [
            f"--{type}={label}=",
            f"import/neo4j_headers/{header_filename}",
            f"import/{target_subdirectory}/*.{target_extension}",
        ]

    def _escape_arg_glob(self):
        for i in self.node_relationship_argument_list:
            i[2] = i[2].replace("*.", ".*\\.")

    def arg_builder(self, input_sg_modules, hmmlist):
        """Reduce the expected sg_module list based on the input/selected modules

        Args:
            input_sg_modules (list): sg_modules to build Neo4j admin import arguments for
            hmmlist (list): if hmms are in input_sg_modules then filter which HMM arguments to create for Neo4j admin import
        """
        self.node_relationship_argument_list = []
        sg_mod_object = SocialgeneModules()
        header_object = Neo4jImportData()
        # retrieve the nodes and relationships that correspond to the input sg_modules
        reduced_dict = {
            "nodes": sg_mod_object.filter_nodes(input_sg_modules),
            "relationships": sg_mod_object.filter_relationships(input_sg_modules),
        }
        for rd_key, rd_value in reduced_dict.items():
            for sg_mod_key, header_keys in rd_value.items():
                if sg_mod_key == "hmms":
                    keylist = [i for i in header_keys if i in hmmlist]
                    keylist.extend(rd_value["base_hmm"])
                else:
                    keylist = header_keys
                for i in keylist:
                    dict_elem = getattr(header_object, rd_key)[i]
                    label = dict_elem["neo4j_label"]
                    header_filename = dict_elem["header_filename"]
                    target_subdirectory = dict_elem["target_subdirectory"]
                    target_extension = dict_elem["target_extension"]
                    self.node_relationship_argument_list.append(
                        self._single_arg_string_builder(
                            label=label,
                            header_filename=header_filename,
                            target_subdirectory=target_subdirectory,
                            target_extension=target_extension,
                            type=rd_key,
                        )
                    )

    def _check_files(self):
        missing_input_headers = []
        missing_input_data = []
        p = Path(self.neo4j_top_dir)
        print(p)
        for i in self.node_relationship_argument_list:
            if not list(p.glob(i[1])):
                # If the files aren't found, check for gzipped versions and change the name
                if list(p.glob(i[1] + ".gz")):
                    i[1] = i[1] + ".gz"
                else:
                    missing_input_headers.append(i[1])
            if not list(p.glob(i[2])):
                # If the files aren't found, check for gzipped versions and change the name
                if list(p.glob(i[2] + ".gz")):
                    i[2] = i[2] + ".gz"
                else:
                    missing_input_data.append(i[2])
        if missing_input_headers or missing_input_data:
            if missing_input_headers:
                message = f"Couldn't find expected header files for neo4j admin import:\n{missing_input_headers}"
            elif missing_input_data:
                message = f"Couldn't find expected data files for neo4j admin import:\n{missing_input_data}"
            else:
                message = "??"
            raise ValueError(message)
        else:
            log.info("Found all expected input files. Huzzah!")

    def _neo4j_admin_import_args(self):
        """Function used to create a list of all the commands/arguments for a docker/neo4j admin import run

        Returns:
            list: commands/arguments to pass to a subprocess
        """

        pre = [
            "docker",
            "run",
            "--rm",
            "--interactive",
            f"--publish={self.dbms_connector_http_listen_address}:{self.dbms_connector_http_listen_address}",
            f"--publish={self.dbms_connector_bolt_listen_address}:{self.dbms_connector_bolt_listen_address}",
            f"--volume={self.neo4j_top_dir}/data:/var/lib/neo4j/data",
            f"--volume={self.neo4j_top_dir}/import:/var/lib/neo4j/import",
            f"--volume={self.neo4j_top_dir}/plugins:/var/lib/neo4j/plugins",
            f"--volume={self.neo4j_top_dir}/import.report:/var/lib/neo4j/import.report",
            f"--user={self.uid}:{self.gid}",
            f"neo4j:{self.neo4j_version}",
            "neo4j-admin",
            "import",
        ]
        # self.node_relationship_argument_list will be inserted between pre and post
        joined_args = [
            f"{i[0]}{i[1]},{i[2]}" for i in self.node_relationship_argument_list
        ]
        post = [
            '--delimiter="\\t"',
            "--high-io=true",
            f"--processors={self.cpus}",
            "--database=neo4j",
            "--ignore-empty-strings=true",
            "--ignore-extra-columns=true",
            "--skip-bad-relationships",
            "--skip-duplicate-nodes",
            # self.additional_args,
        ]
        return pre + list(set(joined_args)) + post

    @staticmethod
    def _run(docker_args):
        run_subprocess(
            command_list=" ".join(docker_args),
            check=False,
            shell=True,
            capture_output=True,
        )

    def run_neo4j_admin_import(self):
        """Run the created commands/arguments in a separate process"""
        self.create_neo4j_directories(self.neo4j_top_dir)
        self.arg_builder(self.input_sg_modules, self.hmmlist)
        self._check_files()
        self._escape_arg_glob()
        self._run(self._neo4j_admin_import_args())
