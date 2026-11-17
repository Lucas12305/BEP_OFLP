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


def get_feas_points(circles, x_lb, x_ub, y_lb, y_ub) -> float:
    # We create 1000 possible facility locations, we then check for each which are feasible to get points to try to fit the number of facilities.
    poss_fac = []
    x_coords = np.linspace(x_lb,x_ub,100)
    y_coords = np.linspace(y_lb,y_ub,100)
    for x in x_coords:
        for y in y_coords:
            p = [x,y]
            feas = True
            for c in circles:
                if IO.interior_circle(p,c):
                    feas = False
                    break
            if feas:
                poss_fac.append(p)
    return poss_fac


def same_point(p, q, eps : float = 1e-9) -> bool:
    return abs(p[0] - q[0]) < eps and abs(p[1] - q[1]) < eps

def find_n_clique(points, fac_fac_dist, n) -> list:
    '''
    This definition finds an n-clique in a graph G. If it finds one it terminates immediatly and returns it.
    :param G: The graph that is considered.
    :param n: The size of the clique that is being searched for.
    '''
    if n <= 0:
        return []
    if n == 1:
        return [0] if len(points) > 0 else []
    if len(points) < n:
        return []
    
    G = nx.Graph()
    pts = np.array(points)
    d = ((pts[:, None, :] - pts[None, :, :]) ** 2).sum(axis=2)
    ii, jj = np.where(np.triu(d >= fac_fac_dist, k=1))
    G.add_edges_from(zip(ii.tolist(), jj.tolist()))
    G = nx.k_core(G, n - 1)

    for v in G:
        if G.degree(v) < n - 1:
            continue
        else:
            clique = [v]
            neighbours = []
            for w in G.neighbors(v):
                if w > v:
                    neighbours.append(w)
            level_neighbours = [neighbours]
            clique_explored = False
            while not clique_explored:
                if len(clique) == n:
                    clique_points = []
                    for i in clique:
                        clique_points.append(points[i])
                    return clique_points
                if len(clique) + len(level_neighbours[-1]) < n:
                    clique.pop()
                    level_neighbours.pop()
                    if len(clique) == 0:
                        clique_explored = True
                else:
                    w = level_neighbours[-1].pop()
                    if G.degree(w) < n - 1:
                        continue
                    poss_vert = []
                    for wi in level_neighbours[-1]:
                        if G.has_edge(w, wi):
                            poss_vert.append(wi)
                    clique.append(w)
                    level_neighbours.append(poss_vert)
    return []

def create_init_sol(facilities, communities, x_lb, x_ub, y_lb, y_ub, fac_fac_dist, r_guess : float=0) -> float:
    # The idea of the algorithm is to decrease the radius until we find a feasible radius, this way we create a good intial solution.
    feas_sol_found = False
    max_iter = 1e3 * len(communities)
    t = 0
    if r_guess == 0:
        r = 2.5 / math.sqrt(len(communities))
    else:
        r = r_guess
    while (not feas_sol_found) and (t < max_iter):
        circles = []
        for i in communities:
            circles.append([i[0], i[1], r])
        points = get_feas_points(circles, x_lb, x_ub, y_lb, y_ub)
        # We now try to fit the facilities in the points while maintaining the fac_fac_distance
        clique = find_n_clique(points, fac_fac_dist, facilities)
        if clique == []:
            r *= 0.95
            t += 1
        else:
            feas_sol_found = True
            return r, clique
    return r, clique

def creating_and_guessing(facilities, communities, x_lb, x_ub, y_lb, y_ub, fac_fac_dist) -> float:
    first_r = create_init_sol(facilities, communities, x_lb, x_ub, y_lb, y_ub, fac_fac_dist)
    r_guess = first_r[0] * 1.01
    i = 0
    n = 20
    while i < n:
        new_r = create_init_sol(facilities, communities, x_lb, x_ub, y_lb, y_ub, fac_fac_dist, r_guess * 1.01)
        if new_r[0] > r_guess:
            r_guess = new_r[0]
            i += 1
        else:
            return r_guess
    return r_guess


def combinations(facilities : int, poss_places : int) -> int:
    '''
    This definition creates a list of lists of all of the possible combinations of facilities and hulls that can be
    :param facilities: The number of facilities.
    :param poss_places: The number of possible places for these facilities.
    '''
    if poss_places <= 0:
        return []

    all_counts = [[]]
    for h in range(poss_places - 1):
        new_counts = []
        for counts in all_counts:
            used = sum(counts)
            for k in range(facilities - used + 1):
                new_counts.append(counts + [k])
        all_counts = new_counts

    result = []
    for counts in all_counts:
        result.append(counts + [facilities - sum(counts)])
    return result



def preliminary_cuts(circles_belong_points : list, circles : list, x_lb : float, x_ub : float, y_lb : float, y_ub : float) -> list:
    '''
    This definition takes a set of feasible points in the form of [x,y, circle, circle/box bound].
    It then sees for every circle the intersections and connects the points that should be connected.
    Lastly, it adds the box extremes to the graph that is the closest.
    :param circles_belong_points: This is a list with for each circle or bound a list with which points it has created.
    :param circles: These are the 4 box bounds and then all of the circles.
    '''
    m = len(circles_belong_points)
    for i in range(4, m):
            print(f"We have a circle {circles[i - 4]} and it has the following points {circles_belong_points[i]}")
    connections = []
    for sp in range(4):
        connections.append([])
        set_of_points = circles_belong_points[sp]
        n = len(set_of_points)
        if n <= 1:
            continue
        else:
            set_of_points.sort(key=lambda  point : (point[0], point[1]))
        for i in range(n - 1):
            p1, p2 = set_of_points[i], set_of_points[i + 1]
            p_middle = [(p2[0] + p1[0]) / 2, (p2[1] + p1[1]) / 2]
            p_middle_inside = False
            for c in circles:
                if IO.interior_circle(p_middle, c):
                    p_middle_inside = True
                    break
            if not p_middle_inside:
                connections[sp].append([p1,p2])


    for sp in range(4,m):
        set_of_points = circles_belong_points[sp]
        connections.append([])
        n = len(set_of_points)
        if n <= 1:
            continue
        else:
            a_p = []
            for p in set_of_points:
                a_p.append([p[0], p[1], IO.angle_points(circles[sp - 4], p1)[0]])
            
            
            # How do we sort, we want to sort such that the angle increases.
            a_p.sort(key=lambda point: point[2])
            for p in range(n):
                p1, p2 = a_p[p], a_p[(p + 1) % n]
                p_middle = IO.angle_to_points(circles[sp - 4], ((p2[2] - p1[2]) / 2 ) % (2 * math.pi))
                # We now check if p_middle is inside some circle
                p_middle_inside = False
                if (x_lb <= p_middle[0] <= x_ub) and (y_lb <= p_middle[1] <= y_ub):
                    for c in circles:
                        if IO.interior_circle(p_middle,c):
                            p_middle_inside = True
                            break
                    if not p_middle_inside:
                        connections[sp].append([p1,p2])

    # We now have a list of all of the sorted connections inside the circles. 
    # This means that we can pick a point and then go both ways, one way is that we go up in index, and the other we go down.
    # This way we handle cases in which a point is made by two circles that intersect exaclty at a point
    subsets_of_hulls = []
    connections_left = []
    for con in connections:
        for c in con:
            connections_left.append(c)

    while len(connections_left) > 0:
        # We take a connection out of the list, so every iteration the list gets shorter
        p1, p2 = connections_left.pop(0)
        set_i = [p1, p2]
        stop_first = IO.at_bound(p1, x_lb, x_ub, y_lb, y_ub)
        stop_last = IO.at_bound(p2, x_lb, x_ub, y_lb, y_ub)

        found = True
        while found and not (stop_first and stop_last):
            found = False
            for con in connections_left[:]:
                if (not stop_last) and same_point(con[0], set_i[-1]):
                    set_i.append(con[1])
                    connections_left.remove(con)
                    found = True
                elif (not stop_first) and same_point(con[1], set_i[0]):
                    set_i.insert(0, con[0])
                    connections_left.remove(con)
                    found = True
            if IO.at_bound(set_i[0], x_lb, x_ub, y_lb, y_ub):
                stop_first = True
            if IO.at_bound(set_i[-1], x_lb, x_ub, y_lb, y_ub):
                stop_last = True
            if same_point(set_i[0], set_i[-1]):
                stop_first = True
                stop_last = True
        subsets_of_hulls.append(set_i)
    print(f"We have derived {len(subsets_of_hulls)} subsets of hulls.")
    # A plot of all of the subsets of hulls, with the number of subsets in the title
    fig, ax = plt.subplots()
    box = patches.Rectangle((x_lb, y_lb), x_ub - x_lb, y_ub - y_lb, linewidth=3, edgecolor='black', facecolor='none')
    ax.add_patch(box)
    for c in circles:
        circle = patches.Circle((c[0], c[1]), c[2], linewidth=1, edgecolor='grey', facecolor='lightgrey', alpha=0.4)
        ax.add_patch(circle)

    colors = plt.cm.tab20.colors
    for k in range(len(subsets_of_hulls)):
        xs = []
        ys = []
        for p in subsets_of_hulls[k]:
            xs.append(p[0])
            ys.append(p[1])
        color = colors[k % len(colors)]
        ax.plot(xs, ys, '-o', color=color, linewidth=2, markersize=4)
        # The number of the subset is written at the middle point of the subset
        mid = subsets_of_hulls[k][len(subsets_of_hulls[k]) // 2]
        ax.text(mid[0], mid[1], str(k), color=color, fontsize=12, fontweight='bold')

    ax.set_xlim(x_lb - abs(x_lb - x_ub) / 10, x_ub + abs(x_lb - x_ub) / 10)
    ax.set_ylim(y_lb - abs(y_lb - y_ub) / 10, y_ub + abs(y_lb - y_ub) / 10)
    ax.set_aspect('equal')
    ax.set_title(f"{len(subsets_of_hulls)} subsets of hulls")
    plt.grid(True)
    plt.show()
    return subsets_of_hulls