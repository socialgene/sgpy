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
        level="NOTSET",
        #  format="%(filename)s/%(module)s/%(funcName)s\::: %(message)s",
        format="%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            RichHandler(rich_tracebacks=True, tracebacks_word_wrap=False, console=c)
        ],
    )

    log = logging.getLogger(__name__)

except ImportError:
    log = logging.getLogger(__name__)

log.setLevel(env_vars["SOCIALGENE_LOGLEVEL"])
# handler = logging.StreamHandler(stream=sys.stdout)

# log.addHandler(handler)
logging.getLogger("neo4j").setLevel(logging.WARNING)
logging.getLogger("asyncio").setLevel(logging.WARNING)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    log.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))


sys.excepthook = handle_exception

# https://blog.hay-kot.dev/fastapi-and-rich-tracebacks-in-development/
# https://coralogix.com/blog/python-logging-best-practices-tips/
