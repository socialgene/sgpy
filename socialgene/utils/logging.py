import logging
import sys

from socialgene.config import env_vars

# If {rich} is installed use it, otherwise.... don't
try:
    from rich.console import Console
    from rich.logging import RichHandler

    c = Console(width=150)
    # https://rich.readthedocs.io/en/stable/logging.html
    logging.basicConfig(
        level=env_vars["SOCIALGENE_LOGLEVEL"],
        #  format="%(filename)s/%(module)s/%(funcName)s\::: %(message)s",
        format="%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            RichHandler(rich_tracebacks=True, tracebacks_word_wrap=False, console=c)
        ],
    )

    log = logging.getLogger("rich")

except ImportError:
    log = logging.getLogger()


logger = logging.getLogger(__name__)
handler = logging.StreamHandler(stream=sys.stdout)

logger.addHandler(handler)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logger.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))


sys.excepthook = handle_exception

# https://blog.hay-kot.dev/fastapi-and-rich-tracebacks-in-development/
# https://coralogix.com/blog/python-logging-best-practices-tips/
