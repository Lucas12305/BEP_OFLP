import re
import sys
import math
import time
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pyscipopt import Model, quicksum
import networkx as nx
from tqdm import tqdm
import math

import Convex_hull as IO


# ## Reading the instances

try:
    BASE_DIR = Path(__file__).parent
except NameError:
    BASE_DIR = Path.cwd()
INSTANCE_DIR = BASE_DIR / "initialized_instances"

SET_RE    = re.compile(r"(\w+)\s+\w+\s*/\s*(\d+)\*(\d+)\s*/")
PARAM_RE  = re.compile(r"community_([xy])\('(\d+)'\)\s*=\s*([-+0-9.eE]+)\s*;")
SCALAR_RE = re.compile(r"Parameter\s+(\w+)\s*/\s*([-+0-9.eE]+)\s*/")
BOUND_RE  = re.compile(r"\b(x|y|facility_facility_sq_distance)\.(lo|up)\([^)]*\)\s*=\s*([^;]+);")


def read_instance(filename: str) -> dict:
    '''Reads a GAMS instance: facilities, communities, box bounds and facility-facility bounds.'''
    text = (INSTANCE_DIR / filename).read_text()

    sets = {name: int(end) - int(start) + 1 for name, start, end in SET_RE.findall(text)}

    coords = {"x": {}, "y": {}}
    for axis, idx, value in PARAM_RE.findall(text):
        coords[axis][int(idx)] = float(value)
    communities = [(coords["x"][k], coords["y"][k]) for k in range(1, sets["j"] + 1)]

    scalars = {name: float(value) for name, value in SCALAR_RE.findall(text)}
    bounds = {}
    for var, side, value in BOUND_RE.findall(text):
        value = value.strip()
        bounds.setdefault(var, {})[side] = scalars.get(value) if value in scalars else float(value)

    return {
        "facilities": sets["i"],
        "x_bounds": (bounds["x"]["lo"], bounds["x"]["up"]),
        "y_bounds": (bounds["y"]["lo"], bounds["y"]["up"]),
        "ff_sq_bounds": (bounds["facility_facility_sq_distance"]["lo"],
                         bounds["facility_facility_sq_distance"]["up"]),
        "communities": communities
    }


def natural_key(path):
    '''Sorts "instance2" before "instance10" by comparing numbers as numbers.'''
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", path.name)]


instances = sorted(
    (f for f in INSTANCE_DIR.iterdir() if f.is_file()),
    key=natural_key,
)


def find_instance(name : str) -> Path:
    '''
    Looks up an instance in the instance folder. The name can be given with or without .gms.
    :param name: The name of the instance, for example ofl_2_100i.
    '''
    for f in instances:
        if name in (f.name, f.stem):
            return f
    raise SystemExit(f"Instance '{name}' was not found in {INSTANCE_DIR}.")


def get_solution(OFP, instance : dict) -> tuple | None:
    '''
    This reads the facilities and the objective L from the best solution of the model.
    It returns None if SCIP did not find a solution.
    :param OFP: The model that has been optimized.
    :param instance: The instance that has been solved, as given by read_instance.
    '''
    if OFP.getNSols() == 0:
        return None
    sol = OFP.getBestSol()
    var = {v.name: v for v in OFP.getVars()}
    facilities = [[OFP.getSolVal(sol, var[f"x_{i}"]), OFP.getSolVal(sol, var[f"y_{i}"])] for i in range(instance['facilities'])]
    L = OFP.getSolVal(sol, var['L'])
    return facilities, L