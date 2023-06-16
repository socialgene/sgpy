import sys

import pkg_resources

from socialgene.utils.logging import log


def find_entrypoints():
    entrypoints = [
        ep.name
        for ep in pkg_resources.iter_entry_points("console_scripts")
        if ep.module_name.startswith("socialgene") and ep.name != "socialgene"
    ]
    return entrypoints


def main():
    log.info("Socialgene has the following entry points for direct command line use:\n")
    log.info("\n".join(find_entrypoints()))


if __name__ == "__main__":
    main()
    sys.exit(main() or 0)
