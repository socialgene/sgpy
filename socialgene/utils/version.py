from importlib.metadata import version

from socialgene.config import env_vars


def main():
    print(f"{version('socialgene')}")


if __name__ == "__main__":
    main()
