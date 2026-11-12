#!/usr/bin/env python
# coding: utf-8

# # Running the obnoxious facility location problem
# This document contains the code to run an instance from the command line, see run.sh.
# The definitions are the ones from Running.ipynb, together with the plot of the solution and a main that reads the command line.

# ## Initial imports

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

import Implementation_OFLP as IO


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


# ## Creating the model

def points_to_cons(points : list, p_inside : list) -> tuple[list, list]:
    A = []
    b = []
    n = len(points)
    # We want constriants of the form Ax <= b, where x = [x,y]
    for i in range(n):
        x1,y1,x2,y2 = points[i][0], points[i][1], points[(i + 1) % n][0], points[(i + 1) % n][1]
        new_A = np.array([y2 - y1 , x1 - x2])
        new_b = np.array(x1 * y2 - x2 * y1)
        if new_A @ p_inside > new_b:
            A.append(-new_A)
            b.append(- new_b)
        else:
            A.append(new_A)
            b.append(new_b)
    return A, b


def feas_region_constraint(x_lb, x_ub, y_lb, y_ub, circles, method : int, model, var) -> None:
    '''
    This function creates the constraint for the feasible region.
    :param method: This is either 0 or 1 and it denotes the methods, 0 is LP, and 1 is GS.
    '''
    if method == 0:
        solutions = IO.generate_points(x_lb, x_ub, y_lb, y_ub)
        b,A = IO.feas_region_box_LP(x_lb, x_ub, y_lb, y_ub, circles, solutions)
    else:
        sorted_points, deleted_points = IO.feas_region_box_GS(x_lb, x_ub, y_lb, y_ub, circles)
        if len(deleted_points) == 0:
            n = len(sorted_points)
            p1, p2 = sorted_points[0], sorted_points[math.ceil(n/2)]
            x1, x2 = p1[0], p1[1]
            y1, y2 = p2[0], p2[1]
            deleted_points = [[x1 + ((x2 - x1) / 2), y1 + ((y2 - y1) / 2)]]
        A, b = points_to_cons(sorted_points, deleted_points[0])
    for i in range(len(b)):
        model.addCons((A[i][0] * var[0] + A[i][1] * var[1]) <= b[i], f'Convex hull constraint {i} for {var}')
    return None


def distance(p1, p2) -> float:
    return (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2


def create_objective(model, var_x, var_y, circles):
    L = model.addVar(vtype='C', lb=0, ub=None, name='L')
    for j, c in enumerate(circles):
        for i, (x, y) in enumerate(zip(var_x, var_y)):
            model.addCons(distance(c,(x,y)) >= L, f"MaxDistance_{i}_With{c[0],c[1]}")
    model.setObjective(L,'maximize')
    return None


def create_LP(x_lb : float, x_ub : float, y_lb : float, y_ub : float, circles : list, facilities : int, min_fac_dist : float, mode : int) -> None:

    OFP = Model(f'OFP problem with {facilities} facilities (with minimal distance {min_fac_dist}) and {len(circles)} communities')
    ObnFacilities_x = [OFP.addVar(vtype='C', lb=x_lb, ub=x_ub, name=f"x_{i}") for i in range(facilities)]
    ObnFacilities_y = [OFP.addVar(vtype='C', lb=y_lb, ub=y_ub, name=f"y_{i}") for i in range(facilities)]
    for i in range(facilities):
        var = [ObnFacilities_x[i], ObnFacilities_y[i]]
        # The minimal distance between the facilities has to be added to OFP, and a facility does not need a distance to itself
        for j in range(i + 1, facilities):
            var2 = [ObnFacilities_x[j], ObnFacilities_y[j]]
            OFP.addCons(distance(var, var2) >= min_fac_dist, f'MinDistF{i}_F{j}')
    create_objective(OFP, ObnFacilities_x,ObnFacilities_y, circles)
    return OFP

def generate_heuristic_solution(OFP, time_limit)-> float:
    '''
    This function creates a first heuristic solution for the instance of the obnoxious facility location problem.
    '''
    OFP.setParam('limits/time', time_limit)
    OFP.optimize()
    sol = OFP.getBestSol()
    Init_sol = OFP.getSolVal(sol, var['L'])
    return Init_sol


# ## Reading and plotting the solution

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


def plot_solution(OFP, instance : dict, save_path : str = None) -> None:
    '''
    This makes a plot of the solution that SCIP has found for the obnoxious facility location problem.
    Around every community it draws a grey circle with the radius of the objective, the facilities are placed outside of the grey area.
    :param OFP: The model that has been optimized.
    :param instance: The instance that has been solved, as given by read_instance.
    :param save_path: If this is given the plot is saved to this file instead of shown.
    '''
    solution = get_solution(OFP, instance)
    if solution is None:
        print("There is no solution to plot.")
        return None
    facilities, L = solution

    x_lb, x_ub = instance['x_bounds']
    y_lb, y_ub = instance['y_bounds']
    communities = instance['communities']

    # Create figure and axes
    fig, ax = plt.subplots()

    # Create a geometric box (Rectangle)
    box = patches.Rectangle((x_lb, y_lb), x_ub - x_lb, y_ub - y_lb, linewidth=3, edgecolor='black', facecolor='none')
    ax.add_patch(box)

    # Plot the communities with a grey circle around them, the facilities have to be outside of the grey area
    # L is a squared distance, so the radius is the square root of L
    # The fill is see-through, so overlapping circles get darker and every circle stays visible
    for c in communities:
        circle = patches.Circle((c[0], c[1]), np.sqrt(L), linewidth=0.8, edgecolor='gray', facecolor=(0.5, 0.5, 0.5, 0.15))
        ax.add_patch(circle)
        ax.scatter(c[0], c[1], s=10, c='black', marker='o', zorder=20)

    # Plot the facilities
    for i, f in enumerate(facilities):
        ax.scatter(f[0], f[1], s=70, c='red', zorder=21)
        ax.annotate(f"F{i}", (f[0], f[1]), textcoords="offset points", xytext=(6, 6), zorder=22)

    # Set limits and equal aspect ratio
    ax.set_xlim(x_lb - abs(x_lb - x_ub) / 2, x_ub + abs(x_lb - x_ub) / 2)
    ax.set_ylim(y_lb - abs(y_lb - y_ub) / 2, y_ub + abs(y_lb - y_ub) / 2)
    ax.set_aspect('equal')
    ax.set_title(f"L = {L:.4f}, radius = {np.sqrt(L):.4f}")

    # Show the plot, or save it when a file is given
    plt.grid(True)
    if save_path is None:
        plt.show()
    else:
        fig.savefig(save_path, dpi=150)
        plt.close(fig)
        print(f"The plot has been saved to {save_path}")
    return None

def update_LP(OFP, circles, x_lb, x_ub, y_lb, y_ub, mode, facilities, time_limit) -> None:
    '''
    This function adds the communities constraints to the model.
    '''
    for var in facilities:
        feas_region_constraint(x_lb, x_ub, y_lb, y_ub, circles, mode, OFP, var)
    return OFP

# ## Running

def run(instance : dict, mode: int, plot : bool = False, plot_file : str = None, presolve_limit : float = 1000, time_limit : float = None) -> None:
    '''
    This functions runs our program
    :param instance: The instance that we want to solve, as given by read_instance.
    :param mode: This is either 0 or 1 and it denotes the methods, 0 is LP, and 1 is GS.
    :param plot: This is a bool which turns on or off the plot of the solution.
    :param plot_file: If this is given the plot is saved to this file instead of shown.
    :param time_limit: The time limit for SCIP in seconds, when it is None there is no limit.
    '''
    facilities = instance['facilities']
    x_lb, x_ub = instance['x_bounds']
    y_lb, y_ub = instance['y_bounds']
    min_FF_dist = instance['ff_sq_bounds'][0]
    communities = instance['communities']
    OFP = create_LP(x_lb, x_ub, y_lb, y_ub, communities, facilities, min_FF_dist, mode)
    Init_sol = generate_heuristic_solution(OFP, presolve_limit)

    circles = []
    for i in communities:
        circles.append([i[0], i[1], Init_sol])
    
    OFP = update_LP(OFP, circles, x_lb, x_ub, y_lb, y_ub, mode, facilities)
    start = time.perf_counter()
    if time_limit is not None:
        OFP.setParam('limits/time', time_limit)
    else:
        OFP.setParam('limits/time', (24 * 3600))
    OFP.optimize()
    OFP.printStatistics()
    OFP.printBestSol()
    total_time = time.perf_counter() - start

    
    # A short summary at the end, so that the result is easy to find in the text file
    print("=" * 60)
    print("Summary")
    print(f"Facilities: {facilities}, communities: {len(communities)}, mode: {mode}")
    print(f"Status: {OFP.getStatus()}, total time: {total_time:.2f} s")
    solution = get_solution(OFP, instance)
    if solution is None:
        print("No solution was found.")
    else:
        facility_points, L = solution
        print(f"Objective L (squared distance): {L}")
        print(f"Radius sqrt(L): {np.sqrt(L)}")
        for i, f in enumerate(facility_points):
            print(f"Facility {i}: x = {f[0]}, y = {f[1]}")
    print("=" * 60, flush=True)

    if plot == True:
        plot_solution(OFP, instance, plot_file)
    return None


def main() -> None:
    '''
    This reads the command line and runs the instance, it is used by run.sh.
    '''
    parser = argparse.ArgumentParser(description="Run an instance of the obnoxious facility location problem.")
    parser.add_argument("instance", help="The name of the instance, for example ofl_2_100i")
    parser.add_argument("--mode", type=int, default=1, choices=[0, 1], help="0 is LP, 1 is GS (default 1)")
    parser.add_argument("--plot-file", default=None, help="Saves the plot of the solution to this file")
    parser.add_argument("--time-limit", type=float, default=None, help="The time limit for SCIP in seconds")
    args = parser.parse_args()

    # When the plot is saved to a file we do not need a window, so this also works without a screen
    if args.plot_file is not None:
        plt.switch_backend("Agg")

    path = find_instance(args.instance)
    instance = read_instance(path.name)
    print(f"Instance: {path.name}", flush=True)
    run(instance, args.mode, plot=args.plot_file is not None, plot_file=args.plot_file, presolve_limit=600,  time_limit=args.time_limit)
    return None


if __name__ == "__main__":
    main()
