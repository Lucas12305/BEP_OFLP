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
                    return clique
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
    with tqdm(total=max_iter, desc="Shrinking r") as pbar:
        while (not feas_sol_found) and (t < max_iter):
            circles = []
            for i in communities:
                circles.append([i[0], i[1], r])
            points = get_feas_points(circles, x_lb, x_ub, y_lb, y_ub)
            # We now try to fit the facilities in the points while maintaining the fac_fac_distance
            clique = find_n_clique(points, fac_fac_dist, facilities)
            if clique == []:
                r *= 0.95
                pbar.update(1)
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
        new_r = create_init_sol(facilities, communities, x_lb, x_ub, y_lb, y_ub, fac_fac_dist, r_guess)
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
    connections = [[]] * m
    for sp in range(1):
        set_of_points = circles_belong_points[sp]
        n = len(set_of_points)
        if n <= 1:
            continue
        else:
            set_of_points.sort(lambda  point : point[1])
        for i in range(n):
            p1, p2 = set_of_points[i], set_of_points[(i + 1) % n]
            p_middle = (p2 + p1) / 2
            p_middle_inside = False
            for c in circles:
                if interior_circle(p_middle, c):
                    p_middle_inside = True
                    break
            if not p_middle_inside:
                connections[i].append([p1,p2])


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
            a_p.sort(key=lambda point: point[2])
            for p in range(n):
                p1, p2 = a_p[p], a_p[(p + 1) % n]
                p_middle = angle_to_points(circles[sp], ((p2[2] - p1[2]) / 2 ) % (2 * math.pi))
                # We now check if p_middle is inside some circle
                p_middle_inside = False
                for c in circles:
                    if interior_circle(p_middle,c):
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
        connections_left.append(con)
    
    while len(connections_left) > 0:
        # In each iteration we pick a random connection [c1, c2]
        cons_con = connections_left[0]
        p1, p2 = cons_con

        # We first check if either are at a bound, then for this ci it stops.
        if at_bound(p1, x_lb, x_ub, y_lb, y_ub):
            stop_iter_clockwise = True
        else:
            stop_iter_clockwise = False
        if at_bound(p2, x_lb, x_ub, y_lb, y_ub):
            stop_iter_anti_clock = True
        else:
            stop_iter_anti_clock = False

        # We look for a connection with c2 as the start and for a connection with c1 as the end. 
        n = 2 * len(connections_left)
        t = 0
        set_i = [p1, p2]
        while (not stop_iter_anti_clock) and (not stop_iter_clockwise) and (t <= n):
            first, last = set_i[0], set_i[-1]
            it = 0
            while ((not found_f) or (not found_l)) and it <= n / 2:
                con = connections_left[it]
                if con[0] == last:
                    set_i.append(con[1])
                    connections_left.remove(con)
                    found_l = True
                elif con[1] == first:
                    set_i.insert(0,con[0])
                    connections_left.remove(con)
                    found_f = True
                it += 1

            # We keep this going until both are at a boundary or if they come to the same point, i.e. connect to make a cycle.
            if at_bound(set_i[0], x_lb, x_ub, y_lb, y_ub):
                stop_iter_clockwise = True
            if at_bound(set_i[-1], x_lb, x_ub, y_lb, y_ub):
                stop_iter_anti_clock = True
            if set_i[0] == set_i[-1]:
                stop_iter_anti_clock = True
                stop_iter_clockwise = True

            # We then found a set, add this to the 
            subsets_of_hulls.append(set_i)
            t += 1



    return subsets_of_hulls