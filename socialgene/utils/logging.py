#!/usr/bin/env python

# python dependencies
import logging
import sys

# external dependencies
# from rich.logging import RichHandler # try/catch below

# internal dependencies


# If {rich} is installed use it, otherwise.... don't
try:
    from rich.logging import RichHandler

    # https://rich.readthedocs.io/en/stable/logging.html
    logging.basicConfig(
        level="NOTSET",
        #  format="%(filename)s/%(module)s/%(funcName)s\::: %(message)s",
        format="%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[RichHandler(rich_tracebacks=True)],
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
