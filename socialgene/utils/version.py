from importlib.metadata import version

from socialgene.config import env_vars


def main():
    print(f"{version('socialgene')}")


def neo4j():
    print(env_vars["NEO4J_VERSION"])


if __name__ == "__main__":
    main()
