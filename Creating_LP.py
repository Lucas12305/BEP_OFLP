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
import Plotting as PLOT
import Precutting as PC



# ## Creating the model

def points_to_cons(points : list, p_inside : list = None) -> tuple[list, list]:
    A = []
    b = []
    n = len(points)
    # We want constriants of the form Ax <= b, where x = [x,y]
    if n < 3:
        return A,b
    if p_inside is None:
        p_inside = [(points[0][0] + points[math.ceil(n/2)][0]) / 2, (p[0][1] + p[math.ceil(n/2)][1]) / 2 ]
    
    for i in range(n):
        x1,y1,x2,y2 = points[i][0], points[i][1], points[(i + 1) % n][0], points[(i + 1) % n][1]
        a1, a2, b1 = y2 - y1 , x1 - x2, x1 * y2 - x2 * y1
        length = math.sqrt(a1 ** 2 + a2 ** 2)
        if length < 1e-12:
            continue
        a1, a2, b1 = a1/length, a2/length, b/length
        if a1 * p_inside[0] + a2 * p_inside[1] > b1:
            a1, a2, b1 = -a1, -a2, -b1
        A.append([a1,a2])
        b.append(b1)
    return A, b


def feas_region_constraint(x_lb, x_ub, y_lb, y_ub, circles, method : int, model, var) -> None:
    '''
    This function creates the constraint for the feasible region when no precutting is used.
    :param method: This is either 0 or 1 and it denotes the methods, 0 is LP, and 1 is GS.
    '''
    if method == 0:
        solutions = IO.generate_points(x_lb, x_ub, y_lb, y_ub)
        planes = IO.feas_region_box_LP(x_lb, x_ub, y_lb, y_ub, circles, solutions)[1]
        A = [[]]*len(planes)
        b = [[]]*len(planes)
        for i in range(len(planes)):
            A[i] = planes[i][:2]
            b[i] = planes[i][2]
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
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)


def create_objective(model, var_x, var_y, circles, best_sol : int=0):
    L = model.addVar(vtype='C', lb=0, ub=None, name='L')
    for j, c in enumerate(circles):
        for i, (x, y) in enumerate(zip(var_x, var_y)):
            model.addCons(((c[0] - x) ** 2 + (c[1] - y) ** 2) >= L, f"MaxDistance_{i}_With{c[0],c[1]}")
    model.setObjective(L,'maximize')
    if best_sol != 0:
        model.addCons(L >= best_sol, f"Using earlier found solution to cut off subproblems")
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
            OFP.addCons(((var[0] - var2[0]) ** 2 + (var[1] - var2[1]) ** 2) >= min_fac_dist, f'MinDistF{i}_F{j}')
    create_objective(OFP, ObnFacilities_x,ObnFacilities_y, circles)
    return OFP


def add_conv_hull_cons_precut_GS(model, var_x, var_y, points) -> None:
    n = len(points)
    p1, p2 = points[0], points[math.ceil(n/2)]
    x1, x2 = p1[0], p1[1]
    y1, y2 = p2[0], p2[1]
    p_ins = [x1 + ((x2 - x1) / 2), y1 + ((y2 - y1) / 2)]
    A,b = points_to_cons(points, p_ins)
    for i in range(len(b)):
        model.addCons((A[i][0] * var_x + A[i][1] * var_y) <= b[i], f'Convex hull constraint {i} for {[var_x,var_y]}')
    return model

def add_conv_hull_cons_precut_LP(model, var_x, var_y, planes) -> None:
    for i in range(len(planes)):
        model.addCons((planes[i][0] * var_x + planes[i][1] * var_y) <= planes[i][2], f'Convex hull constraint {i} for {[var_x,var_y]}')
    return model

def update_LP(OFP, circles, x_lb, x_ub, y_lb, y_ub, mode, facilities) -> None:
    '''
    This function adds the communities constraints to the model.
    '''
    # The variables are stored as x_0..x_(n-1), y_0..y_(n-1), L, so facility i has x at index i and y at index n + i
    all_vars = OFP.getVars()
    for i in range(facilities):
        feas_region_constraint(x_lb, x_ub, y_lb, y_ub, circles, mode, OFP, [all_vars[i], all_vars[facilities + i]])
    return OFP



