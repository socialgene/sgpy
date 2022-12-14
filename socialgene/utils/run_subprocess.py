# python dependencies
import subprocess

# external dependencies
from rich.console import Console

# internal dependencies
from socialgene.utils.logging import log

console = Console()


def run_subprocess(
    command_list,
    check=True,
    shell=False,
    **kwargs,
):
    """Run something in a separate process

    Args:
        command_list (list or str): list if shell=False, str if shell=True # TODO:better preventative
        check (bool, optional): [description]. Defaults to True.
        shell (bool, optional): [description]. Defaults to False.

    Raises:
        ValueError: [description]
    """
    if not command_list:
        log.error("Attempted to pass an empty command to subprocess.run()")
        raise ValueError
    if isinstance(command_list, list):
        command_list_string = " ".join(command_list)
    else:
        command_list_string = command_list
    log.info(f"Executing external program:\n{command_list_string}")
    with console.status(
        "",
        spinner="bouncingBar",
    ) as status:
        subprocess.run(command_list, check=check, shell=shell, **kwargs)
        _ = status
