#!/usr/bin/env python
# coding: utf-8

# # Running the obnoxious facility location problem
# This document contains the code to run an instance from the command line, see run.sh.
# The definitions are the ones from Running.ipynb, together with the plot of the solution and a main that reads the command line.
import re
import sys
import math
import time
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pyscipopt import Model, Sepa, SCIP_RESULT
import networkx as nx
from tqdm import tqdm
import math
import Convex_hull as IO
import Plotting as PLOT
import Precutting as PC
import Creating_LP as CLP
import Administration_and_loading as AL
import Seperator as SEP

# ## Running

def add_seperator(model, var_x : list, var_y : list, circles : list, root_hulls : list, mode : int, seperator : int, update : int, hull_of_facility : list = None) -> None:
    '''
    This function checks if the seperator has to be added or not
    '''
    if seperator == 0 and update == 0:
        return None
    sepa = SEP.HullSeperator(var_x, var_y, circles, root_hulls, seperator, update, hull_of_facility, mode)
    model.includeSepa(sepa, "hull_sepa", "Constraints of the convex hull of the feasible region", priority=1000, freq=1)
    return None

def print_summary(OFP, instance : dict, total_time : float, settings : str) -> None:
    '''
    A short summary at the end, so that the result is easy to find in the text file.
    '''
    print("=" * 60)
    print("Summary")
    print(f"Facilities: {instance['facilities']}, communities: {len(instance['communities'])}, {settings}")
    print(f"Status: {OFP.getStatus()}, total time: {total_time:.2f} s")
    solution = AL.get_solution(OFP, instance)
    if solution is None:
        print("No solution was found.")
    else:
        facility_points, L = solution
        print(f"Objective L (squared distance): {L}")
        print(f"Radius sqrt(L): {np.sqrt(L)}")
        for i, f in enumerate(facility_points):
            print(f"Facility {i}: x = {f[0]}, y = {f[1]}")
    print("=" * 60, flush=True)
    return None


def run(instance : dict, mode: int, plot : bool = False, plot_file : str = None, presolve_limit : float = 1000, time_limit : float = None, precutting : int = 0, seperator : int = 0, update : int = 0) -> None:
    '''
    This functions runs our program
    :param instance: The instance that we want to solve, as given by read_instance.
    :param mode: This is either 0 or 1 and it denotes the methods to make the hull, 0 is LP, and 1 is GS.
    :param plot: This is a bool which turns on or off the plot of the solution.
    :param plot_file: If this is given the plot is saved to this file instead of shown.
    :param time_limit: The time limit for SCIP in seconds, when it is None there is no limit.
    :param precutting: This is either 0 or 1, 0 is off and 1 is on.
    :param seperator: This is either 0 or 1, 0 is off and 1 is on. When it is on the hull only gets added as constraints that the LP solution violates.
    :param update: 0 does not update the hull during the process, 1 makes the hull again from scratch (LP or GS) with the new bounds, 2 cuts the hull of the last node (only GS).
    '''
    assert update in (0, 1) or mode == 1, "Updating from the last hull is only possible with GS (mode 1)."
    facilities = instance['facilities']
    x_lb, x_ub = instance['x_bounds']
    y_lb, y_ub = instance['y_bounds']
    min_FF_dist = instance['ff_sq_bounds'][0]
    communities = instance['communities']
    Init_sol = PC.create_init_sol(facilities, communities, x_lb, x_ub, y_lb, y_ub, min_FF_dist)
    circles = []
    for i in communities:
        circles.append([i[0], i[1], Init_sol[0]])
    settings = f"mode: {mode}, precutting: {precutting}, seperator: {seperator}, update: {update}"
    solutions = IO.generate_points(x_lb, x_ub, y_lb, y_ub)
    start = time.perf_counter()


    if precutting == 0:
        # This means that we only have one box that we have shrunken using the LP method or the GS method
        OFP = CLP.create_LP(x_lb, x_ub, y_lb, y_ub, communities, facilities, min_FF_dist, mode)
        if seperator == 0 and update == 0:
            OFP = CLP.update_LP(OFP, circles, x_lb, x_ub, y_lb, y_ub, mode, facilities)
        else:
            # The constraints of the hull are not in the model, the seperator adds them
            if mode == 0:
                root_hull = IO.feas_region_box_LP(x_lb, x_ub, y_lb, y_ub, circles, solutions)[1]
            else:
                root_hull = list(IO.feas_region_box_GS(x_lb, x_ub, y_lb, y_ub, circles))
            all_vars = OFP.getVars()
            add_seperator(OFP, all_vars[:facilities], all_vars[facilities:2 * facilities], circles, [root_hull], mode, seperator, update)
        OFP.setParam('limits/time', 3600 if time_limit is None else time_limit)
        OFP.optimize()
        OFP.printStatistics()
        OFP.printBestSol()
        print_summary(OFP, instance, time.perf_counter() - start, settings)
        if plot == True:
            PLOT.plot_solution(OFP, instance, plot_file)
        return None


    else: # Precutting == 1
        # We divide into the different sets
        if mode == 0:
            subsets_of_points, planes = IO.feas_region_box_LP_cutting(x_lb, x_ub, y_lb, y_ub, circles, solutions)
            root_hulls = planes
        else:
            subsets_of_points, not_used_points = IO.feas_region_box_GS_cutting(x_lb, x_ub, y_lb, y_ub, circles)
            root_hulls = [[subsets_of_points[h], not_used_points[h]] for h in range(len(subsets_of_points))]
        poss_places = len(subsets_of_points)
        subproblems_combinations = PC.combinations(facilities, poss_places, subsets_of_points, min_FF_dist, x_lb, x_ub, y_lb, y_ub)
        best_found_sol = Init_sol[0] ** 2
        index_bfs = [-1]

        # We now want to create all of the subproblem
        print(f"We are going to be looking at {len(subproblems_combinations)} combinations. After the best possible solution will be pasted.")
        OFP = {}
        OFP_optimal = {}
        pbar = tqdm(range(len(subproblems_combinations)), desc="Subproblems")
        for i in pbar:
            comb = subproblems_combinations[i]
            hull_of_facility = []
            OFP[i] = Model(f'OFP_{i} problem with {facilities} facilities (with minimal distance {min_FF_dist}) and {len(circles)} communities and combinations {comb}.')
            ObnFacilities_x = [OFP[i].addVar(vtype='C', lb=x_lb, ub=x_ub, name=f"x_{j}") for j in range(facilities)]
            ObnFacilities_y = [OFP[i].addVar(vtype='C', lb=y_lb, ub=y_ub, name=f"y_{j}") for j in range(facilities)]
            for l in range(facilities):
                var = [ObnFacilities_x[l], ObnFacilities_y[l]]
                # The minimal distance between the facilities has to be added to OFP, and a facility does not need a distance to itself
                for k in range(l + 1, facilities):
                    var2 = [ObnFacilities_x[k], ObnFacilities_y[k]]
                    OFP[i].addCons(((var[0] - var2[0]) ** 2 + (var[1] - var2[1]) ** 2) >= min_FF_dist, f'MinDistF{l}_F{k}')

                # Now comes the extra step, here we specifically make it so that the facilities are in a certain spot
                for m in range(poss_places):
                    if comb[m] != 0:
                        comb[m] -= 1
                        hull_of_facility.append(m)
                        if mode == 0:
                            CLP.add_conv_hull_cons_precut_LP(OFP[i], var[0], var[1], root_hulls[m])
                        else:
                            CLP.add_conv_hull_cons_precut_GS(OFP[i], var[0], var[1], subsets_of_points[m])
                        break

            # The hulls stay in the model, the seperator makes them smaller when the bounds of a node get smaller
            add_seperator(OFP[i], ObnFacilities_x, ObnFacilities_y, circles, root_hulls, mode, seperator, update, hull_of_facility)
            CLP.create_objective(OFP[i], ObnFacilities_x, ObnFacilities_y, circles, best_found_sol)
            OFP[i].setParam('limits/time', 3600 if time_limit is None else time_limit)
            OFP[i].hideOutput()
            OFP[i].optimize()
            if OFP[i].getNSols() > 0:
                L = OFP[i].getObjVal()
                if L > best_found_sol:
                    best_found_sol = L
                    index_bfs = [i]
                    for OFPI in OFP_optimal.values():
                        OFPI.freeProb()
                    OFP_optimal = {i : OFP[i]}
                elif abs(L - best_found_sol) < 1e-9:
                    index_bfs = index_bfs + [i]
                    OFP_optimal[i] = OFP[i]
                else:
                    OFP[i].freeProb()
            else:
                OFP[i].freeProb()

            del OFP[i]
            pbar.set_postfix(best_radius=f"{np.sqrt(best_found_sol):.4f}")

        for j in index_bfs:
            if j == -1:
                print(f"The best solution found was the initials solution: {Init_sol}, no subproblem found a better one.")
                continue
            OFP_optimal[j].printStatistics()
            OFP_optimal[j].printBestSol()

        if index_bfs[0] == -1:
            print("=" * 60)
            print("Summary")
            print(f"Facilities: {facilities}, communities: {len(communities)}, {settings}")
            print(f"Total time: {time.perf_counter() - start:.2f} s")
            print("No subproblem found a solution better than the initial solution.")
            print(f"Objective L (squared distance): {best_found_sol}")
            print(f"Radius sqrt(L): {np.sqrt(best_found_sol)}")
            print("=" * 60, flush=True)
            return None

        print_summary(OFP_optimal[index_bfs[0]], instance, time.perf_counter() - start, settings)
        if plot == True:
            PLOT.plot_solution(OFP_optimal[index_bfs[0]], instance, plot_file)
        return None
            



def main() -> None:
    '''
    This reads the command line and runs the instance, it is used by run.sh.
    '''
    parser = argparse.ArgumentParser(description="Run an instance of the obnoxious facility location problem.")
    parser.add_argument("instance", help="The name of the instance, for example ofl_2_100i")
    parser.add_argument("--mode", type=int, default=1, choices=[0, 1], help="0 is LP, 1 is GS (default 1)")
    parser.add_argument("--precutting", type=int, default=0, choices=[0, 1], help="0 is off, 1 is on (default 0)")
    parser.add_argument("--seperator", "--separator", dest="seperator", type=int, default=0, choices=[0, 1], help="0 is off, 1 is on (default 0)")
    parser.add_argument("--update", type=int, default=0, choices=[0, 1, 2], help="0 does not update the hull, 1 makes it again from scratch, 2 cuts the hull of the last node, only with GS (default 0)")
    parser.add_argument("--plot-file", default=None, help="Saves the plot of the solution to this file")
    parser.add_argument("--time-limit", type=float, default=None, help="The time limit for SCIP in seconds")
    args = parser.parse_args()

    # When the plot is saved to a file we do not need a window, so this also works without a screen
    if args.plot_file is not None:
        plt.switch_backend("Agg")

    path = AL.find_instance(args.instance)
    instance = AL.read_instance(path.name)
    print(f"Instance: {path.name}", flush=True)
    run(instance, args.mode, plot=args.plot_file is not None, plot_file=args.plot_file, presolve_limit=600, time_limit=args.time_limit, precutting=args.precutting, seperator=args.seperator, update=args.update)
    return None


if __name__ == "__main__":
    main()



