import sys
from importlib.metadata import entry_points


def find_entrypoints():
    entrypoints = [
        i.name
        for i in entry_points(group="console_scripts")
        if i.value.startswith("socialgene")
    ]
    entrypoints.sort()
    return entrypoints


def main():
    print("Socialgene has the following entry points for direct command line use:")
    print("\n".join(find_entrypoints()))


if __name__ == "__main__":
    main()
    sys.exit(main() or 0)
