import subprocess

from socialgene.utils.logging import CONSOLE, log


def run_subprocess(
    command_list,
    check=True,
    shell=False,
    capture_output=True,
    status=False,
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
    if status:
        log.info(f"Executing external program:\n{command_list_string}")
        with CONSOLE.status(
            "",
            spinner="bouncingBar",
        ) as status:
            result = subprocess.run(
                command_list,
                check=check,
                shell=shell,
                capture_output=capture_output,
                **kwargs,
            )
    else:
        log.debug(f"Executing external program:\n{command_list_string}")
        result = subprocess.run(
            command_list,
            check=check,
            shell=shell,
            capture_output=capture_output,
            **kwargs,
        )
    if result.returncode != 0:
        log.error(f"Error code: {result.returncode}")
        log.error(f"Error stdout: {result.stderr.decode('utf-8')}")
        raise SystemExit
