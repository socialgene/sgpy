#!/usr/bin/env python


import os

import pkg_resources

os_vars = dict(os.environ)
env_vars = {}

env_file = pkg_resources.resource_filename(__name__, "common_parameters.env")
internal_vars = {}
with open(env_file) as f:
    for line in f:
        if line.startswith("#") or line == "\n":
            continue
        k, v = line.strip().split("=", 1)
        internal_vars[k] = v


for k, v in internal_vars.items():
    if isinstance(v, str):
        try:
            internal_vars[k] = int(v)
        except ValueError:
            internal_vars[k] = v


for k, v in internal_vars.items():
    if isinstance(v, str):
        try:
            internal_vars[k] = float(v)
        except ValueError:
            internal_vars[k] = v


for k, v in internal_vars.items():
    if v == "NONE":
        internal_vars[k] = None


for k, v in internal_vars.items():
    if k in os_vars:
        env_vars[k] = os_vars[k]
    else:
        env_vars[k] = v
