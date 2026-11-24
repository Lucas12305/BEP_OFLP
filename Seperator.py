from pyscipopt import Sepa, SCIP_RESULT
import Convex_hull as IO
import Creating_LP as CLP

def box_key(x_lb, x_ub, y_lb, y_ub, digits=9):
    return (round(x_lb, digits), round(x_ub, digits), round(y_lb, digits), round(y_ub, digits))

def facet_key(a : float,b : float,c : float, rounding : int=9):
    return (round(a,rounding), round(b,rounding), round(c,rounding))


class HullSeperator(Sepa):
    def __init__(self, var_x, var_y, circles, root_hulls, sep_mode=1, update=0, hull_deriving_mode=0, hull_of_facility=None):
        self.var_x, self.var_y, self.circles = var_x, var_y, circles, A, b
        self.sep_mode, self.update, self.hull_deriving_mode = sep_mode, update, hull_deriving_mode
        self.root_hulls = root_hulls
        if hull_of_facility is None:
            hull_of_facility = [0] * len(var_x)
        self.hull_of_facility = hull_of_facility


        self.root_facets = []
        for points, not_used in root_hulls:
            A, b = CLP.points_to_cons(points)
            self.root_facets.append({facet_key(A[i][0], A[i][1], b[i]) for i in range(len(b))})

        # The hulls that we stored per box and the box that belongs to each of them
        self.hulls = {}
        self.hull_boxes = {}
        self.var_index = None
        # The constraints that are already added, so that we do not add them again in the next round of the same node
        self.done = set()
        self.root_done = set()


    def cut_hull(self, hull : list, xy_cut : int, cut_val : float, bound_type : int) -> list:
        '''
        Cuts a hull once. bound_type is 0 when cut_val is a new lower bound and 1 when it is a new upper bound.
        '''
        points, not_used = hull
        if len(points) < 3:
            return [points, not_used]

        all_p = points + not_used
        lo = min(p[xy_cut] for p in all_p)
        hi = max(p[xy_cut] for p in all_p)

        # The cut only does something when it lies inside the hull
        if bound_type == 0 and not cut_val > lo + 1e-9:
            return [points, not_used]
        if bound_type == 1 and not cut_val < hi - 1e-9:
            return [points, not_used]
        # When the cut lies outside of all of the points there is nothing feasible left
        if (bound_type == 0 and cut_val > hi) or (bound_type == 1 and cut_val < lo):
            return [[], []]

        other = 1 - xy_cut
        box = [0, 0, 0, 0]
        box[2 * xy_cut], box[2 * xy_cut + 1] = lo, hi
        box[2 * other] = min(p[other] for p in all_p)
        box[2 * other + 1] = max(p[other] for p in all_p)

        b1, not_used1, b2, not_used2 = IO.cutting_convex_hull(box[0], box[1], box[2], box[3], points, not_used, self.circles, xy_cut, cut_val)
        if bound_type == 0:
            return [b2, not_used2]
        return [b1, not_used1]


    def get_cutval(self, node) -> list:
        branchings = node.getParentBranchings()
        if branchings is None:
            return None
        variables, cut_vals, bound_types = branchings
        var, cut_val, bound_type = variables[-1], cut_vals[-1], bound_types[-1]

        # We look up which facility and which coordinate the variable is, using the index of the variable
        if self.var_index is None:
            self.var_index = {}
            for j in range(len(self.var_x)):
                self.var_index[self.model.getTransformedVar(self.var_x[j]).getIndex()] = (j, 0)
                self.var_index[self.model.getTransformedVar(self.var_y[j]).getIndex()] = (j, 1)
        if var.getIndex() not in self.var_index:
            return None
        facility, xy_cut = self.var_index[var.getIndex()]
        return facility, xy_cut, cut_val, bound_type
   
        
    def get_hull(self, x_lb : float, x_ub : float, y_lb : float, y_ub : float, hull_deriving_mode : int, update : int, h : int):
        '''
        :param hull_deriving_mode: 0 is that we don't update hull, 1 is making hull from scratch, 2 is updating from earlier
        '''
        bounds = (x_lb, x_ub, y_lb, y_ub)
        root = self.root_hulls[h]
        if hull_deriving_mode == 0:
            return root

        if hull_deriving_mode == 1:
            hull = self.cut_hull(root, bounds)
        
        else: # Mode == 2
            key = (h,) + box_key(x_lb, x_ub, y_lb, y_ub)
            if key in self.hulls:
                self.n_stored_hits += 1
                return self.hulls[key]
            # A hull that we stored for bounds that contain these bounds, contains all of the feasible points of these bounds.
            # We take the one with the fewest points, because that is the one that is cut the most.
            start_hull = root
            for old_key, old_bounds in self.hull_boxes.items():
                if old_key[0] != h:
                    continue
                if old_bounds[0] <= x_lb + 1e-9 and old_bounds[1] >= x_ub - 1e-9 and old_bounds[2] <= y_lb + 1e-9 and old_bounds[3] >= y_ub - 1e-9:
                    cand = self.hulls[old_key]
                    if len(cand[0]) + len(cand[1]) < len(start_hull[0]) + len(start_hull[1]):
                        start_hull = cand
            hull = self.cut_hull(start_hull, bounds)
            self.hulls[key] = hull
            self.hull_boxes[key] = bounds

        return hull

    def add_row(self, x, y, a1 : float, a2 : float, c : float, name : str, local : bool, force : bool) -> bool:
        '''
        Adds the constraint a1 * x + a2 * y <= c to the LP, and gives True when this makes the node infeasible.
        '''
        model = self.model
        row = model.createEmptyRowSepa(self, name=name, lhs=None, rhs=c, local=local)
        model.cacheRowExtensions(row)
        model.addVarToRow(row, x, a1)
        model.addVarToRow(row, y, a2)
        model.flushRowExtensions(row)
        is_feas = model.addCut(row, forcecut=force)
        model.releaseRow(row)
        return is_feas

    def sepaexeclp(self):
        model = self.model
        node = model.getCurrentNode()
        found = False 
        for j in range(len(self.var_x)):
            x = model.getTransformedVar(self.var_x[j])
            y = model.getTransformedVar(self.var_y[j])

            
            h = self.hull_of_facility[j] # starting hull of facility
            box = (x.getLbLocal(), x.getUbLocal(), y.getLbLocal(), y.getUbLocal())
            points, not_used = self.get_hull(box[0], box[1], box[2], box[3], self.hull_deriving_mode, self.update, h)
            A, b = CLP.points_to_cons(points)

        
            # When mode == 1 we use the seperator, otherwise the seperator gets used to update the hull after a cut
            if self.sep_mode == 1:
                px, py = x.getLPSol(), y.getLPSol()

                # We track how much each constraint is violated
                violated = []
                for i in range(len(A)):
                    violation = A[i][0] * px + A[i][1] * py - b[i]
                    if violation > 1e-6:
                        violated.append((violation, i))

                # We now pick the 5 most violated ones
                violated.sort(reverse=True)
                for violation, i in violated[:5]:
                    if self.add_row(x,y,A[i][0], A[i][1], b[i], f"hull_f{j}_c{i}", True, False):
                        return {"result": SCIP_RESULT.CUTOFF}
                    found = True
            
            else: #self.sep_mode == 0
                # We add each constraint that is new
                for i in range(len(A)):
                    key = facet_key(A[i][0], A[i][1], b[i])
                    if key in self.root_facets[h]:
                        if node.getDepth() == 0 and j not in self.root_done:
                            if self.add_row(x, y, A[i][0], A[i][1], b[i], f"hull_{j}_{i}", False, True):
                                return {"result": SCIP_RESULT.CUTOFF}
                            found = True
                    
                    # The new constraints only hold in this node and below, and we do not add them twice in one node
                    elif (node.getNumber(), j, key) not in self.done:
                        self.done.add((node.getNumber(), j, key))
                        if self.add_row(x, y, A[i][0], A[i][1], b[i], f"hull_f{j}_c{i}", True, True):
                            return {"result": SCIP_RESULT.CUTOFF}
                        found = True
                if node.getDepth() == 0:
                    self.root_done.add(j)
                    
        if found:
            return {"result": SCIP_RESULT.SEPARATED}
        return {"result": SCIP_RESULT.DIDNOTFIND}             

