# python dependencies
import os
import pwd
from pathlib import Path

# external dependencies

# internal dependencies
from socialgene.config import env_vars
from socialgene.neo4j.schema.socialgene_modules import SocialgeneModules
from socialgene.utils.run_subprocess import run_subprocess
from socialgene.utils.logging import log


class Neo4jAdminImport(SocialgeneModules):
    """Socialgene's Neo4jAdminImport class handles building the Docker/Neo4j admin import command line arguments"""

    def __init__(
        self,
        neo4j_top_dir: str,
        cpus: int = 1,
        additional_args_to_pass_to_neo4j: str = "",
        uid: int = None,
        gid: int = None,
        module_list: list = None,
        hmm_list: list = None,
        dbms_connector_http_listen_address: int = 7474,
        dbms_connector_bolt_listen_address: int = 7687,
        neo4j_version=env_vars["NEO4J_VERSION"],
        *args,
        **kwargs,
    ):
        """Class for building and executing Neo4j admin import

        Args:
            neo4j_top_dir (str): directory containing the necessary dirs and files for Neo4j admin import
            module_list (list): name of grouped nodes/relationships (socialgene/neo4j/schema/define_modules.py)
            hmm_list (list): list of hmms source databases used (socialgene/neo4j/schema/define_hmmlist.py)
            cpus (int, optional): [description]. Defaults to 1.
            additional_args_to_pass_to_neo4j (str, optional): [description]. Defaults to "".
            uid (int, optional): [description]. Defaults to None.
            gid (int, optional): [description]. Defaults to None.
        """
        super().__init__(*args, **kwargs)
        self.input_sg_modules = module_list
        self.input_hmmlist = hmm_list
        self.neo4j_top_dir = neo4j_top_dir
        self.cpus = cpus
        self.additional_args_to_pass_to_neo4j = additional_args_to_pass_to_neo4j
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
        # return f"--{type}={label}=import/neo4j_headers/{header_filename},import/{target_subdirectory}/^.*{chr(92)}.{target_extension}.*"
        return [
            f"--{type}={label}=",
            f"import/neo4j_headers/{header_filename}",
            f"import/{target_subdirectory}/*.{target_extension}.*",
        ]

    def _escape_arg_glob(self):
        for i in self.node_relationship_argument_list:
            i[2] = i[2].replace("*.", ".*\\.")

    def get_nodes_and_relationships(self, module_list, hmm_list):
        """Reduce the expected sg_module list based on the input/selected modules

        Args:
            module_list (list): sg_modules to build Neo4j admin import arguments for
            hmm_list (list): if hmms are in module_list then filter which HMM arguments to create for Neo4j admin import
        """
        self.add_modules(module_list)
        self.add_hmms(hmm_list)

    def build_nodes_and_relationships_argument_list(self):
        for node in self.nodes:
            self.node_relationship_argument_list.append(
                self._single_arg_string_builder(
                    label=node.neo4j_label,
                    header_filename=node.header_filename,
                    target_subdirectory=node.target_subdirectory,
                    target_extension=node.target_extension,
                    type="node",
                )
            )
        for rel in self.relationship:
            self.node_relationship_argument_list.append(
                self._single_arg_string_builder(
                    label=rel.neo4j_label,
                    header_filename=rel.header_filename,
                    target_subdirectory=rel.target_subdirectory,
                    target_extension=rel.target_extension,
                    type="node",
                )
            )

    def _check_files(self):
        missing_input_headers = []
        missing_input_data = []
        p = Path(self.neo4j_top_dir)
        print(p)
        for i in self.node_relationship_argument_list:
            if not list(p.glob(i[1])):
                missing_input_headers.append(i[1])
            if not p.glob(i[2]):
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

    def _neo4j_admin_import_args(self, docker=False):
        """Function used to create a list of all the commands/arguments for a docker/neo4j admin import run

        Returns:
            list: commands/arguments to pass to a subprocess
        """

        docker_command = [
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
        ]

        neo4j_admin_command = ["neo4j-admin", "database", "import", "full"]

        # self.node_relationship_argument_list will be inserted between pre and post
        joined_args = [
            f"{i[0]}{i[1]},{i[2]}" for i in self.node_relationship_argument_list
        ]
        post = [
            '--delimiter="\\t"',
            "--high-parallel-io=On",
            f"--threads={self.cpus}",
            "--ignore-empty-strings=true",
            "--ignore-extra-columns=true",
            "--skip-bad-relationships",
            "--skip-duplicate-nodes",
            "neo4j",  # this is the database name. If I remember correctly this must be "neo4j" for the community edition of neo4j
        ]

        if docker:
            pre = docker_command + neo4j_admin_command
        else:
            pre = neo4j_admin_command
        return pre + list(set(joined_args)) + post

    @staticmethod
    def _run(cmd):
        run_subprocess(
            command_list=" ".join(cmd),
            check=False,
            shell=True,
            capture_output=True,
        )

    def run_neo4j_admin_import(self, docker=False):
        """Run the created commands/arguments in a separate process"""
        self.create_neo4j_directories(self.neo4j_top_dir)
        self.get_nodes_and_relationships(self.input_sg_modules, self.hmmlist)
        self.build_nodes_and_relationships_argument_list()
        self._check_files()
        self._escape_arg_glob()
        self._run(self._neo4j_admin_import_args(docker))
