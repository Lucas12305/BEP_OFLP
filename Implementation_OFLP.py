#!/usr/bin/env python
# coding: utf-8

# # Implementation obnoxious facility location problem
# This document contains the cleaned up code for the addition proposed in the Bachelor End Project of Lucas Klaasse. 

# # Backend definitions

# ## Initial imports
# These packages are used in the definitions that follow.

# In[ ]:


import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import math
from pyscipopt import Model, quicksum


# ## General definitions
# The following definitions are general definitions or implementations of well known formulas

# In[ ]:


def solve_abc_form(a,b,c) -> list:
    '''
    This function is an implementation of the abc formula for solving x in a*x^2 + b*x + c = 0.
    :param a: a
    :param b: b
    :param c: c
    '''
    D = (b ** 2) - (4 * a * c) 
    if D < 0:
        return []
    else:
        return [float((- b + np.sqrt(D)) / (2 * a)), float(( - b - np.sqrt(D)) / (2 * a))]



def interior_circle(point : list, circle : list) -> bool:
    '''
    This function looks whether or not a point is inside the circle.
    :param point: The point that we want to check.
    :param circle: The circle that we want to check.
    '''
    if (circle[0] - point[0]) ** 2 + (circle[1] - point[1]) ** 2 < circle[2] ** 2 - 1e-4:
        return True
    else:
        return False


def generate_points(x_lb : float, x_ub : float, y_lb : float, y_ub : float) -> list:
    '''
    A function that creates points inside the box for testing the implementation of the convex hull creation with the LP.
    :param x_lb: The lowest x of the box.
    :param x_ub: The highest x of the box.
    :param y_lb: The lowest y of the box.
    :param y_ub: The highest y of the box.
    '''
    points = []
    x = np.linspace(x_lb, x_ub, 100)
    y = np.linspace(y_lb, y_ub, 100)
    for i in x:
        for j in y:
            points.append([i,j])
    return points


def angle_points(refP : list, oP : list) -> tuple:
    '''
    Takes a reference point refP and calculates the angel in radians between the refP and the other point oP and the distance.
    :param refP: The reference point.
    :param oP: The other point.
    '''
    dx = oP[0] - refP[0]
    dy = oP[1] - refP[1]    
    radians = float(np.arctan2(dy, dx) % (2 * np.pi))
    distance = float(np.sqrt(dy ** 2 + dx ** 2)) 
    return radians, distance


def angle_to_points(circle : list, angle : float) -> list:


# ## Calculating intersections
# The following definitions are for calculating intersections in the main program of the code

# In[ ]:


def box_circles(x_lb, x_ub, y_lb, y_ub, circles : list, index_c : int, circles_belong_point : list) -> list:
    '''
    This function takes a box and a circle and returns the intersections that this circle has with the box.
    :param x_lb: The lowest x of the box.
    :param x_ub: The highest x of the box.
    :param y_lb: The lowest y of the box.
    :param y_ub: The highest y of the box.
    :param circle: The circle that we are considering
    '''
    x, y, r = circles[index_c]
    intersections = []

    # We first check which sides the circle intersects with
    if x - r <= x_lb <= x + r:
        a = 1
        b = -2 * y
        c = x ** 2 - 2 * x * x_lb + x_lb ** 2 + y ** 2 - r ** 2
        py = solve_abc_form(a,b,c)
        if py != []:
            for yi in py:
                p = [x_lb, yi]
                if (p not in intersections) and y_lb <= yi <= y_ub:
                    intersections.append(p)
                if (p not in circles_belong_point[0]) and y_lb <= yi <= y_ub:
                    circles_belong_point[0].append(p)
    
    if x - r <= x_ub <= x + r:
        a = 1
        b = -2 * y
        c = x ** 2 - 2 * x * x_ub + x_ub ** 2 + y ** 2 - r ** 2
        py = solve_abc_form(a,b,c)
        if py != []:
            for yi in py:
                p = [x_lb, yi]
                if (p not in intersections) and y_lb <= yi <= y_ub:
                    intersections.append(p)
                if (p not in circles_belong_point[1]) and y_lb <= yi <= y_ub:
                    circles_belong_point[1].append(p)
    
    if y - r <= y_lb <= y + r:
        a = 1
        b = -2 * x
        c = y ** 2 - 2 * y * y_lb + y_lb ** 2 + x ** 2 - r ** 2
        px = solve_abc_form(a,b,c)
        if px != []:
            for xi in px:
                p = [x_lb, yi]
                if (p not in intersections) and x_lb <= xi <= x_ub:
                    intersections.append(p)
                if (p not in circles_belong_point[2]) and x_lb <= xi <= x_ub:
                    circles_belong_point[2].append(p)
            
    if y - r <= y_ub <= y + r:
        a = 1
        b = -2 * x
        c = y ** 2 - 2 * y * y_ub + y_ub ** 2 + x ** 2 - r ** 2
        px = solve_abc_form(a,b,c)
        if px != []:
            for xi in px:
                p = [x_lb, yi]
                if (p not in intersections) and x_lb <= xi <= x_ub:
                    intersections.append(p)
                if (p not in circles_belong_point[3]) and x_lb <= xi <= x_ub:
                    circles_belong_point[3].append(p)
    
    return intersections, circles_belong_point





def circles_circles(circles : list, c1_index : int, c2_index : int) -> list:
    '''
    This function takes two circles and returns there intersection points.
    :param c1: The first circle considered.
    :param c2: The second circle considered.
    '''
    x0, y0, r0 = circles[c1_index]
    x1, y1, r1 = circles[c2_index]

    dx = x1 - x0
    dy = y1 - y0
    dist = np.sqrt(dx**2 + dy**2)

    # No intersection
    if dist > r0 + r1:
        return []

    # One circle completely inside the other
    if dist < abs(r0 - r1):
        return []

    # Same centre
    if dist == 0:
        return []

    # Distance from centre 0 to midpoint between intersections
    a = (r0**2 - r1**2 + dist**2) / (2 * dist)

    h_squared = r0**2 - a**2
    h_squared = max(0, h_squared)
    h = np.sqrt(h_squared)

    # Midpoint
    xm = x0 + a * dx / dist
    ym = y0 + a * dy / dist

    # Offset from midpoint
    rx = -dy * h / dist
    ry = dx * h / dist

    p1 = [float(xm + rx), float(ym + ry)]
    p2 = [float(xm - rx), float(ym - ry)]
    
    # Tangent circles: only one intersection
    if abs(h) < 1e-4:
        return [p1]

    return [p1, p2]




# ## LP definitions
# These definitions are there to create the convex hull of the set of points using the LP.

# In[ ]:

def update_arcs(arcs : list, TwoPoints : list) -> list:
    '''
    This function takes arcs and new points and then connects the points that seem feasible.
    '''
    if len(arcs) == 0:
        return [TwoPoints[1], TwoPoints[0]]
    for arc in arcs:
        if (arc[0] <= TwoPoints[0]) and (TwoPoints[1] <= arc[1]):
            return None
    
    return None

        

def preliminary_cuts(feasible_points, circles_belong_points : list, circles : list) -> list[list]:
    '''
    This definition takes a set of feasible points in the form of [x,y, circle, circle/box bound].
    It then sees for every circle the intersections and connects the points that should be connected.
    Lastly, it adds the box extremes to the graph that is the closest.
    :param feasible_points: All of the feasible points and where they come from.
    :param circles_belong_points: This is a list with for each circle or bound a list with which points it has created.
    :param circles: These are the 4 box bounds and then all of the circles.
    '''
    # We have for each side and for each circle points that belong to the circle.
    # We want to see which points we have to connect and which we don't
    # We do this with the following
    m = len(circles_belong_points)
    connections = [[]]*m
    for sp in range(4,m):
        set_of_points = circles_belong_points[sp]
        n = len(set_of_points)
        if n <= 1:
            continue
        else:
            a_p = []
            for p in set_of_points:
                a_p.append([p[0], p[1], angle_points(circles[sp], p1)[0]])
            
            
            # How do we sort, we want to sort such that the angle increases.
            a_p.sort()
            for p in range(n):
                p1, p2 = a_p[p], a_p[(p + 1) % n]
                p_middle = angle_to_points(( (p2[2] - p1[2]) / 2 ) % (2 * math.pi), circles[sp])
                # We now check if p_middle is inside some circle
                for c in circles:
                    if not interior_circle(p,c)
                        p_middle_inside = True
                        break
                if not p_middle_inside:
                    connections[sp].append([p1,p2])

    # We now have a list of all of the connections inside the circles. 
    # We now want to connect the points such that we create seperate shapes.
    sets = [feasible_points[0]]]
    for p in range(1, len(feasible_points)):
        numb_sets = len(sets)
        add_to_set = None
        for set in sets[numb_sets]:
            if p in set:
                add_to_set = set
        if add_to_set == None:
            sets.append(p)
            add_to_set = sets[-1]

        # We have 2 cases
        # Case 1: p already belongs to a set, we then search for an addition to this set because of this point p.
        # Subcase 1a : There is one extra point, we add it
        # Subcase 1b : There are 2 or more points, this means that in that point some circles coincide in such a way that more points are generated
        #              We copy the set in the amount of sets with points there are with one point added to each one
        # Case 2: p does not belong to a set so then we add it and try to find 2 points to each side of the point.







   

        

    
        
        




    return Sets




def create_conv_region(points: list, solutions: list) -> list:
    '''
    For each solution, find a cutting plane that separates
    the solution from all points.
    :param points: The points that need to be on one side of the cutting plane.
    :param solutions: The solutions of the program that we want to try to cut.
    '''
    cutting_planes = []
    t = 0

    for s in solutions:
        model = Model("CuttingPlane")
        model.hideOutput()

        # Variables
        p0 = model.addVar("p0", lb=-1, ub=1)
        p1 = model.addVar("p1", lb=-1, ub=1)
        p2 = model.addVar("p2", lb=-1, ub=1)

        # Objective: maximize p0*s[0] + p1*s[1] + p2
        model.setObjective(
            p0 * s[0] + p1 * s[1] + p2,
            sense="maximize"
        )

        # Constraints: p0*x + p1*y + p2 <= 0
        for x, y in points:
            model.addCons(p0*x + p1*y + p2 <= 0)

        # Solve
        model.optimize()

        if model.getStatus() == "optimal":
            plane = [
                model.getVal(p0),
                model.getVal(p1),
                model.getVal(p2)
            ]

            if model.getObjVal() >= 1e-4 and plane not in cutting_planes:
                cutting_planes.append(plane)
                print(f"{s} at t = {t}, {plane}")
                t += 1

    return cutting_planes


# ## GS definitions
# These definitions are there to create the convex hull of the set of points using GS.

# In[ ]:


def eval_clockwise(p1, p2, p3) -> bool:
    '''
    Takes 3 points and evaluates if the line between the second and third point 
    is a clockwise turn from the line between 1st and 2nd point using the cross product.
    :param p1: The first point.
    :param p2: The second point.
    :param p3: The third point.
    '''
    area = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0])
    return area > 0



def create_conv_region_GS(points : list) -> list:
    '''
    An implementation of the Graham scan algorithm to create the convex hull of a set of points.
    :points param: This is a list of points in the (x,y).
    :circles param: These are all of the circles that have made the points.
    '''
    
    # We first find the point with the lowest y-value and set it as the start of the sorted points.
    # When two points have the same y-value the one with the lowest x-value gets taken
    points_GS = points.copy()
    if len(points_GS) == 0:
        return [], []
    if len(points_GS) <= 2:
        return points_GS, []
    min_y = min(points_GS, key=lambda p: (p[1], p[0]))
    sorted_points = [min_y]
    points_GS.remove(min_y)

    # We now calculate the angles between each 2 points
    for p in points_GS:
        angle, dist = angle_points(min_y, p)
        p.append(angle)
        p.append(dist)
    
    points_GS.sort(key=lambda p: (p[2], p[3]))
    # We add the first point and delete it from points
    sorted_points.append(points_GS[0])
    points_GS.remove(points_GS[0])
    del_points = []

    for c in points_GS:
        while len(sorted_points) >= 2:
            a = sorted_points[-2]
            b = sorted_points[-1]
    
            if eval_clockwise(a, b, c):
                break
            else:
                del_points.append(sorted_points.pop()[:2])
    
        sorted_points.append(c)
        
    sorted_points = [p[:2] for p in sorted_points]
    return sorted_points, del_points


def angle_points_custom(p1, p2, direction):
    '''
    This function calculates the angle between 2 points and then turns it so that the right starting point gets used.
    :param p1: The first point considered, i.e. the reference point.
    :param p2: The second point considered, i.e. the other point.
    '''
    
    dx = p2[0] - p1[0]
    dy = p2[1] - p1[1]

    angle = math.atan2(dy, dx)

    reference_angles = [
        math.pi,            # 0: Up
        math.pi/2,          # 1: Left
        0,                  # 2: Down
        3*math.pi / 2       # 3: right
    ]

    if direction not in range(4):
        raise ValueError("Direction must be 0, 1, 2 or 3")
    distance = float(np.sqrt(dy ** 2 + dx ** 2))
    return (angle - reference_angles[direction]) % (2 * math.pi), distance




def create_conv_region_GS_know_points(known_cycle : list, points : list, rotation : int) -> list:
    '''
    An implementation of the Graham scan algorithm to create the convex hull of a set of points when we already know the order of a set of points.
    :param known_cycle: This is a list of the points we already know are in the right order in (x,y).
    :param points: This is a list of points in the (x,y).
    :param rotation: The side that gets picked as the reference side. 0 for x_lb, 1 for y_ub, 2 for x_ub, 3 for y_lb.
    '''
    # We first find the point with the lowest y-value and set it as the start of the sorted points.
    # When two points have the same y-value the one with the lowest x-value gets taken
    points_GS = []
    if len(points) == 0:
        return known_cycle, []
    
    start_vertex = known_cycle[0]
    end_vertex = known_cycle[-1]

    # Angle of the reference line


    for p in points:
        angle, dist = angle_points_custom(start_vertex, p, rotation)
        points_GS.append([p, angle, dist])

    # Sort by adjusted angle
    points_GS.sort(key=lambda p: (p[1], p[2]))
    points_GS1 = points_GS
    
    sorted_points = known_cycle.copy()
    del_points = []
    
    # Protect all known vertices
    min_length = len(known_cycle)

    for p, angle, dist in points_GS1:
        while len(sorted_points) >= min_length:
            a = sorted_points[-2]
            b = sorted_points[-1]
                    
            if eval_clockwise(a, b, p):
                break

            del_points.append(sorted_points.pop())
    
        sorted_points.append(p)

    # Remove duplicated closing vertex
    for i in range(len(sorted_points) - 1):
        if sorted_points[i + 1] == sorted_points[0]:
            return sorted_points[:i+1], del_points
    
    return sorted_points, del_points


# ## Visualization
# The following definitions are used for visualization of what is going on inside the definitions.

# In[ ]:


def plot_box_LP(x_lb,x_ub,y_lb,y_ub, circles, intersections, planes : list, trc :bool = False, rc : list=[]) -> None:
    '''
    This makes a plot for the box that gets cut using the LP.
    '''
    # Create figure and axes
    fig, ax = plt.subplots()

    # Create a geometric box (Rectangle)
    # Arguments: (lower-left x, y), width, height
    box = patches.Rectangle((x_lb, y_lb), x_ub-x_lb, y_ub - y_lb, linewidth=3, edgecolor='black', facecolor='none')
    ax.add_patch(box)
        
    # Create the circles
    # Arguments: (center x, y), radius
    for c in circles:
        circle = patches.Circle((c[0], c[1]), c[2], linewidth=1, edgecolor='black', fill = False)
        ax.add_patch(circle)
    # Add shapes to the axes

    for p in intersections:
        ax.scatter(p[0], p[1], s=20, c='red', marker='o')
    if trc == True:
        for c in rc:
            circle = patches.Circle((c[0], c[1]), c[2], linewidth=1, edgecolor='lightgray', fill = False)
            ax.add_patch(circle)
            ax.scatter(c[0], c[1], s=1, c='blue', marker='o', label='Center Point')
        
    for plane in planes:
        p0, p1, p2 = plane
    
        if abs(p1) > 1e-10:
            x_vals = np.linspace(x_lb, x_ub, 200)
            y_line = -(p0 * x_vals + p2) / p1
    
            # Keep everything inside the box
            y_line = np.clip(y_line, y_lb, y_ub)
    
            if p1 < 0:
                # Shade below line
                ax.fill_between(
                    x_vals, y_lb, y_line,
                    alpha=1, color = 'blue'
                )
            else:
                # Shade above line
                ax.fill_between(
                    x_vals, y_line, y_ub,
                    alpha=1, color = 'blue'
                )
    
            ax.plot(x_vals, y_line, color = 'blue')
    
        else:
            # Vertical line
            if abs(p0) > 1e-10:
                x_line = -p2 / p0
    
                if x_lb <= x_line <= x_ub:
                    ax.vlines(x_line, y_lb, y_ub, color='blue')
    
                    if p0 < 0:
                        # x <= x_line
                        ax.axvspan(x_lb, x_line, alpha=1, color = 'blue')
                    else:
                        # x >= x_line
                        ax.axvspan(x_line, x_ub, alpha=1, color = 'blue')
        
    # Set limits and equal aspect ratio
    ax.set_xlim(x_lb - abs(x_lb-x_ub)/2, x_ub + abs(x_lb-x_ub)/2)
    ax.set_ylim(y_lb - abs(y_lb-y_ub)/2, y_ub + abs(y_lb-y_ub)/2)
    ax.set_aspect('equal')
        
    # Show the plot
    plt.grid(True)
    plt.show()
    return None



def plot_box_GS(x_lb,x_ub,y_lb,y_ub, circles, intersections, sorted_points : list, trc :bool = False, rc : list=[]) -> None:
    '''
    This makes a plot for the box that gets cut using the LP.
    '''
    # Create figure and axes
    fig, ax = plt.subplots()

    # Create a geometric box (Rectangle)
    # Arguments: (lower-left x, y), width, height
    box = patches.Rectangle((x_lb, y_lb), x_ub-x_lb, y_ub - y_lb, linewidth=3, edgecolor='black', facecolor='none')
    ax.add_patch(box)
        
    # Create the circles
    # Arguments: (center x, y), radius
    for c in circles:
        circle = patches.Circle((c[0], c[1]), c[2], linewidth=1, edgecolor='black', fill = False)
        ax.add_patch(circle)
    
    # Add shapes to the axes
    for p in intersections:
        #print(f"We have a point {p}")
        ax.scatter(p[0], p[1], s=70, c='green', zorder=21)
    
    if trc == True:
        for c in rc:
            circle = patches.Circle((c[0], c[1]), c[2], linewidth=1, edgecolor='lightgray', fill = False)
            ax.add_patch(circle)
            ax.scatter(c[0], c[1], s=1, c='blue', marker='o', label='Center Point')
        
    for i in range(len(sorted_points)):
        p1 = sorted_points[i]
        p2 = sorted_points[(i + 1) % len(sorted_points)]
        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            color='blue'
        )
  
    # Set limits and equal aspect ratio
    ax.set_xlim(x_lb - abs(x_lb-x_ub)/2, x_ub + abs(x_lb-x_ub)/2)
    ax.set_ylim(y_lb - abs(y_lb-y_ub)/2, y_ub + abs(y_lb-y_ub)/2)
    ax.set_aspect('equal')
        
    # Show the plot
    plt.grid(True)
    plt.show()
    return None


# ## Creating a convex hull
# These definitions are made to create a convex hull with either the LP method or the GS method.

# In[ ]:


def feas_region_box_LP(x_lb : float, x_ub : float, y_lb : float, y_ub : float, circles : list[list], solutions : list ,epsilon : float = 1e-4, trc : bool = False, plot : bool = False) -> None:
    '''
    This definition takes a box and circles and returns the smallest convex polygon inside the box.
    :param x_lb: lowest x of the box.
    :param x_ub: highest x of the box.
    :param y_lb: lowest y of the box.
    :param y_ub: highest y of the box.
    :param circles: A list of cirles, which are denoted as lists with (x,y,r)
    :param trc: This is a bool which stand for track removed circles, this is to check if no circles are removed that cut of an additional piece.
    :param plot: This is a bool which turns on or off the plot function.
    '''
    assert x_lb < x_ub or y_lb < y_ub, "Make sure that the upperbounds are higher than the lower bounds."
    for circle in circles:
        assert len(circle) == 3, "Make sure all circles have 3 parameters."    
    rc = []
    
    # The following code first looks which circles are (partly inside our space), it could add some circles which are actually not in there
    # due to it just using the bounds, but this is minimal.
    circles_used = []
    for circle in circles:
        if (circle[0] + circle[2] <= x_lb - epsilon) or (circle[0] - circle[2] >= x_ub + epsilon) or (circle[1] - circle[2] >= y_ub + epsilon) or (circle[1] + circle[2] <= y_lb - epsilon):
            if trc == True:
                rc.append(circle)
            continue
        else:
            if not circle in circles_used:
                circles_used.append(circle)
                

    # The following code checks if there are circles inside circles that we then don't have to consider.
    n = len(circles_used)
    remove_indices = []
    for i in range(n):
        for j in range(i + 1,n):
            c1 = circles_used[i]
            c2 = circles_used[j]
            if (j not in remove_indices) and ((((c1[0] - c2[0]) ** 2) + ((c1[1] - c2[1]) ** 2)) ** (1/2)) + (c2[2] ) <= (c1[2]):
                remove_indices.append(j)
            elif (i not in remove_indices) and ((((c1[0] - c2[0]) ** 2) + ((c1[1] - c2[1]) ** 2)) ** (1/2)) + (c1[2] ) <= (c2[2]):
                remove_indices.append(i)
    remove_indices.sort(reverse=True)
    for index in remove_indices:
        if trc == True:
                rc.append(circles_used[index])
        circles_used.remove(circles_used[index])
    
    # We now look at all of the intersections between the circles and the square
    m = len(circles_used) + 4
    points = [(x_lb, y_lb), (x_lb, y_ub), (x_ub, y_lb), (x_ub, y_ub)]
    points_belong_circle = [[]]*m
    for i in range(len(circles_used)):
        point, points_belong_circle = box_circles(x_lb, x_ub, y_lb, y_ub, circles_used, i, points_belong_circle)
        #print(f"We are looking at the circle: {circle} and it produces points :{point}")
        for p in point:
            #print(f"Point {p} is not in the list of points: {points}? {p not in points}")
            if p not in points_belong_circle[i+4]:
                points_belong_circle[i+4].append(p)
            if p not in points:
                points.append(p)

    # We now look at all of the intersections between the circles
    for i in range(n):
        for j in range(i + 1,n):
            c1, c2 = circles_used[i], circles_used[j]
            if c1 == c2:
                continue
            elif ((c1[0] - c2[0]) ** 2) + ((c1[1] - c2[1]) ** 2) >= c1[2] ** 2 + c2[2] ** 2:
                continue
            else:
                intersections = circles_circles(circles,i,j)
            for point in intersections:
                if point not in points:
                    interior = False
                    if (x_lb <= point[0] <= x_ub) and (y_lb <= point[1] <= y_ub):
                        for c in circles_used:
                            if interior_circle(point, c) == True:
                                interior = True
                        if (not interior):
                            points_belong_circle[i+4].append(point)
                            points_belong_circle[j+4].append(point)
                            points.append(point)


    #We now create the new convex feasible region
    planes = create_conv_region(points, solutions)
    
    if plot == True:
        plot_box_LP(x_lb,x_ub,y_lb,y_ub, circles, points, planes, trc, rc)
    
    #return None
    return points, planes





def feas_region_box_GS(x_lb : float, x_ub : float, y_lb : float, y_ub : float, circles : list[list] ,epsilon : float = 1e-4, trc : bool = False, plot : bool = False) -> None:
    '''
    This definition takes a box and circles and returns the smallest convex polygon inside the box.
    :param x_lb: lowest x of the box.
    :param x_ub: highest x of the box.
    :param y_lb: lowest y of the box.
    :param y_ub: highest y of the box.
    :param circles: A list of cirles, which are denoted as lists with (x,y,r)
    :param trc: This is a bool which stand for track removed circles, this is to check if no circles are removed that cut of an additional piece.
    :param plot: This is a bool which turns on or off the plot function.
    '''
    assert x_lb < x_ub or y_lb < y_ub, "Make sure that the upperbounds are higher than the lower bounds."
    for circle in circles:
        assert len(circle) == 3, "Make sure all circles have 3 parameters."    
    rc = []
    
    # The following code first looks which circles are (partly inside our space), it could add some circles which are actually not in there
    # due to it just using the bounds, but this is minimal.
    circles_used = []
    for circle in circles:
        if (circle[0] + circle[2] <= x_lb - epsilon) or (circle[0] - circle[2] >= x_ub + epsilon) or (circle[1] - circle[2] >= y_ub + epsilon) or (circle[1] + circle[2] <= y_lb - epsilon):
            if trc == True:
                rc.append(circle)
            continue
        else:
            if not circle in circles_used:
                circles_used.append(circle)
                

    # The following code checks if there are circles inside circles that we then don't have to consider.
    n = len(circles_used)
    remove_indices = []
    for i in range(n):
        for j in range(i + 1,n):
            c1 = circles_used[i]
            c2 = circles_used[j]
            if (j not in remove_indices) and ((((c1[0] - c2[0]) ** 2) + ((c1[1] - c2[1]) ** 2)) ** (1/2)) + (c2[2] ) <= (c1[2]):
                remove_indices.append(j)
            elif (i not in remove_indices) and ((((c1[0] - c2[0]) ** 2) + ((c1[1] - c2[1]) ** 2)) ** (1/2)) + (c1[2] ) <= (c2[2]):
                remove_indices.append(i)
    remove_indices.sort(reverse=True)
    for index in remove_indices:
        if trc == True:
                rc.append(circles_used[index])
        circles_used.remove(circles_used[index])
    
    # We now look at all of the intersections between the circles and the square
    m = len(circles_used) + 4
    points = [(x_lb, y_lb), (x_lb, y_ub), (x_ub, y_lb), (x_ub, y_ub)]
    points_belong_circle = [[]]*m
    for i in range(len(circles_used)):
        point, points_belong_circle = box_circles(x_lb, x_ub, y_lb, y_ub, circles_used, i, points_belong_circle)
        #print(f"We are looking at the circle: {circle} and it produces points :{point}")
        for p in point:
            #print(f"Point {p} is not in the list of points: {points}? {p not in points}")
            if p not in points_belong_circle[i+4]:
                points_belong_circle[i+4].append(p)
            if p not in points:
                points.append(p)

    # We now look at all of the intersections between the circles
    for i in range(n):
        for j in range(i + 1,n):
            c1, c2 = circles_used[i], circles_used[j]
            if c1 == c2:
                continue
            elif ((c1[0] - c2[0]) ** 2) + ((c1[1] - c2[1]) ** 2) >= c1[2] ** 2 + c2[2] ** 2:
                continue
            else:
                intersections = circles_circles(circles,i,j)
            for point in intersections:
                if point not in points:
                    interior = False
                    if (x_lb <= point[0] <= x_ub) and (y_lb <= point[1] <= y_ub):
                        for c in circles_used:
                            if interior_circle(point, c) == True:
                                interior = True
                        if (not interior):
                            points_belong_circle[i+4].append(point)
                            points_belong_circle[j+4].append(point)
                            points.append(point)

                        
    #We now create the new convex feasible region
    pch, dp = create_conv_region_GS(points)
    
    if plot == True: 
        plot_box_GS(x_lb,x_ub,y_lb,y_ub, circles, points, pch, trc, rc)
    
    return pch, dp


# ## Smart cuts
# These definitions are for making smart cuts in the convex hulls.

# In[ ]:


def box_making(sorted_points : list, first_vert_index1 : int, last_vert_index1 : int, first_vert_index2 : int, last_vert_index2 : int):
    '''
    This function takes the first and last index for each of the boxes and creates the two parts of the cycles.
    :param sorted_points: The list of the old cycle of points that make up the convex hull.
    :param first_vert_index1: The index of the first vertex of b1.
    :param first_vert_index2: The index of the first vertex of b2.
    :param last_vert_index1: The index of the last vertex of b1.
    :param last_vert_index2: The index of the last vertex of b2.
    '''
    if first_vert_index1 >= last_vert_index1:
        b1 = sorted_points[first_vert_index1:] + sorted_points[:last_vert_index1]
    else:
        b1 = sorted_points[first_vert_index1:first_vert_index2]

    if first_vert_index2 >= last_vert_index2:
        b2 = sorted_points[first_vert_index2:] + sorted_points[:last_vert_index2] 
    else:
        b2 = sorted_points[first_vert_index2:last_vert_index2]

    return b1, b2


def get_boundary_index(sorted_points : list, cut_val : float, xy_cut : int) -> float:
    '''
    This definition sorts the sorted cycle for the convex hull as in left, right, below or above the cut value. It does this using the bisection method.
    :param sorted_points: This is the cycle of the points that make up the convex hull.
    :param cut_val: The value that we have made a cut at.
    :param xy_cut: Indicates whether we cut over x (0) or over y (1).
    '''
    assert len(sorted_points) >= 2, "The list 'sorted points' must contain at least 2 elements."
    not_found = True
    max_iter = len(sorted_points)
    n = 0
    i = max_iter  / 2
    sets_i = [[-1, math.floor(i - 1/2)], [math.ceil(i), max_iter - 1]]
    found_interval = False
    while (not found_interval) and (n <= (max_iter + 10)):
        found = False
        n += 1
        for interval in sets_i:
            if (sorted_points[interval[0]][xy_cut] < cut_val) and (sorted_points[interval[1]][xy_cut] >= cut_val):
                i = (interval[1] - interval[0]) / 2
                if i < 0.9:
                    first_i = interval[1]
                    found_interval = True
                else:
                    sets_i = [[interval[0], math.floor(i - 1/2)], [math.ceil(i), interval[1]]]
                found = True
                break
        if not found:
            split_problem = []
            for interval in sets_i:
                i = (interval[1] - interval[0]) / 2
                split_problem.append([interval[0], math.floor(i - 1/2)])
                split_problem.append([math.ceil(i), interval[1]])
            sets_i = split_problem  



    found_interval = False
    n = 0
    max_iter = len(sorted_points)
    j = max_iter / 2
    sets_j = [[-1, math.floor(j - 1/2)], [math.ceil(j), max_iter - 1]]
    while (not found_interval) and (n <= max_iter + 10):
        found = False
        n += 1
        for interval in sets_j:
            if (sorted_points[interval[0]][xy_cut] > cut_val) and (sorted_points[interval[1]][xy_cut] <= cut_val):
                j = (interval[1] - interval[0]) / 2
                if j < 0.9:
                    first_j = interval[1]
                    found_interval = True
                else:
                    sets_j = [[interval[0], math.floor(j - 1/2)], [math.ceil(j), interval[1]]]
                found = True
                break
        if not found:
            split_problem = []
            for interval in sets_j:
                j = (interval[1] - interval[0]) / 2
                split_problem.append([interval[0], math.floor(j - 1/2)])
                split_problem.append([math.ceil(j), interval[1]])
            sets_j = split_problem   
        
    if first_i[xy_cut] == cut_val:
        last_j = first_i
    else:
        last_j = (first_i - 1) % max_iter
    
    if first_j[xy_cut] == cut_val:
        last_i = first_j
    else:
        last_i = (first_j - 1) % max_iter
    
    return box_making(sorted_points, first_i, last_i, first_j, last_j)






def cutting_convex_hull(x_lb : float, x_ub : float, y_lb : float, y_ub : float, sorted_points : list, not_used_points : list, circles : list, xy_cut : int, cut_val : float, plot : bool=False, in_out : bool = False) -> list:
    '''
    This definition takes a box with circles, a convex hull and a cut and then makes a better convex hull for the 2 new boxes.
    :param bounds: the bounds of the box.
    :param sorted_points: These are the points that are used for the convex hull in order, starting with lowest point. 
    :param circles: The circles that we need to consider.
    :param xy_cut: This parameter indicates if we cut over y or x, x = 0, y = 1.
    :param cut_val: The value at which we want to cut.
    '''
    assert xy_cut == 0 or xy_cut == 1, f"Choose either 0 for a cut over x and 1 for a cut over y, not {xy_cut}."
    if xy_cut == 0:
        assert x_lb <= cut_val <= x_ub, f"Choose a cut value between {x_lb} and {x_ub}"
    if xy_cut == 1:
        assert y_lb <= cut_val <= y_ub, f"Choose a cut value between {y_lb} and {y_ub}"


    lb = [x_lb, y_lb]
    ub = [x_ub, y_ub]
    
    # We then calculate all of the intersection that the cut has with the circles and the upper and lower bounds.
    feas_p = []
    v = abs(xy_cut - 1)
    for c in circles:
        x,y,r = c
        if c[xy_cut] - r - 1e-4 <= cut_val <= c[xy_cut] + r + 1e-4:
            a = 1
            b = -2 * c[v]
            c = c[xy_cut] ** 2 - 2 * c[xy_cut] * cut_val + cut_val ** 2 + c[v] ** 2 - r ** 2
            s = solve_abc_form(a,b,c) 
            for si in s:
                if lb[v] <= si <= ub[v]:
                    if xy_cut == 0:
                        pi = [cut_val, si]
                    else:
                        pi = [si, cut_val]
                    in_circ = True
                    for c in circles:
                        if interior_circle(pi, c):
                            in_circ = False
                    if in_circ and pi not in feas_p:
                        feas_p.append(pi)
    
    # We now check if the intersections of the cut and bounds is feasible and if so we add it
    feas_top, feas_bot = True, True
    m = len(circles)
    t = 0
    if xy_cut == 0:
        p_t = [cut_val, y_ub]
        p_b = [cut_val, y_lb]
    else:
        p_t = [x_ub, cut_val] 
        p_b = [x_lb, cut_val]
    while t < m and (feas_top == True or feas_bot == True):
        if interior_circle(p_t, circles[t]):
            feas_top = False
        if interior_circle(p_b, circles[t]):
            feas_bot = False
        t += 1
    if feas_top:
        feas_p.append(p_t)
    if feas_bot:
        feas_p.append(p_b)

    # We now create b1 and b2 using highest and lowest
    # We also make it so that it starts for b1 at the highest and for b2 at the lowest.
    b1, b2 = get_boundary_index(sorted_points, cut_val, xy_cut)
   

    # We also split the not used points in b1 and b2:
    all_items_with_dubbels = not_used_points + feas_p
    all_items = []
    for i in all_items_with_dubbels:
        if i not in all_items and i not in sorted_points:
            all_items.append(i)

    oP1, oP2 = [], []
    for p in all_items:
        if (p[xy_cut] <= cut_val) and (p not in oP1):
            oP1.append(p)
        if (p[xy_cut] >= cut_val) and (p not in oP2):
            oP2.append(p)
    

    #We have 2 types of rotation, one for the x cut and one for a y cut
    if xy_cut == 0:
        rotation = [0,2]
    else:
        rotation = [3,1]

    if in_out:
            print(f"The input for box 1 is: {b1, oP1, rotation[0]}")
            print(f"The input for box 2 is: {b2, oP2, rotation[1]}")

    b1_cycle, b1_not_used = create_conv_region_GS_know_points(b1, oP1, rotation[0]) 
    b2_cycle, b2_not_used = create_conv_region_GS_know_points(b2, oP2, rotation[1]) 

    b1 = b1_cycle
    b2 = b2_cycle
    if xy_cut == 0:
        b1_x_ub = cut_val
        b1_y_ub = y_ub
        b2_x_lb = cut_val
        b2_y_lb = y_lb

    else:
        b1_x_ub = x_ub
        b1_y_ub = cut_val
        b2_x_lb = x_lb
        b2_y_lb = cut_val
    
    if plot == True:
        plot_box_GS(x_lb,b1_x_ub,y_lb,b1_y_ub, circles, b1_not_used, b1, [], [])
        plot_box_GS(b2_x_lb,x_ub,b2_y_lb,y_ub, circles, b2_not_used, b2, [], [])
    
    if in_out:
        print(f"Box 1 has the following cycle {b1}, and has {b1_not_used}")
        print(f"Box 2 has the following cycle {b2}, and has {b2_not_used}")
    return b1, b1_not_used, b2, b2_not_used

