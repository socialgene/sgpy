#!/usr/bin/env python

# python dependencies
import pkg_resources
import os

# external dependencies

# internal dependencies

env_vars = dict(os.environ)

if "HMMSEARCH_IEVALUE" not in env_vars:
    env_file = pkg_resources.resource_filename(__name__, "common_parameters.env")
    env_vars = {}
    with open(env_file) as f:
        for line in f:
            if line.startswith("#") or line == "\n":
                continue
            key, value = line.strip().split("=", 1)
            env_vars[key] = value

for key, value in env_vars.items():
    if isinstance(value, str):
        try:
            env_vars[key] = int(value)
        except ValueError:
            env_vars[key] = value
    if isinstance(value, str):
        try:
            env_vars[key] = float(value)
        except ValueError:
            env_vars[key] = value
    if value == "NONE":
        env_vars[key] = None
