from pyscipopt import Sepa, SCIP_RESULT
import Convex_hull as IO
import Creating_LP as CLP

def box_key(x_lb, x_ub, y_lb, y_ub, digits=9):
    return (round(x_lb, digits), round(x_ub, digits), round(y_lb, digits), round(y_ub, digits))

def facet_key(a : float,b : float,c : float, rounding : int=9):
    return (round(a,rounding), round(b,rounding), round(c,rounding))


class HullSeperator(Sepa):
    def __init__(self, var_x, var_y, circles, root_hulls, sep_mode=1, hull_deriving_mode=0, hull_of_facility=None):
        self.var_x, self.var_y, self.circles = var_x, var_y, circles
        self.sep_mode, self.hull_deriving_mode = sep_mode, hull_deriving_mode
        self.root_hulls = root_hulls
        if hull_of_facility is None:
            hull_of_facility = [0] * len(var_x)
        self.hull_of_facility = hull_of_facility


        self.root_facets = []
        for points, not_used in root_hulls:
            A, b = CLP.points_to_cons(points)
            self.root_facets.append({facet_key(A[i][0], A[i][1], b[i]) for i in range(len(b))})

        # The hull of every facility in every node that we have seen, with as key (node number, facility)
        self.node_hulls = {}
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


    def get_cutval(self, node) -> tuple:
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
   
        
    def get_hull(self, node, j : int, bounds : tuple, hull_deriving_mode : int, h : int) -> list:
        '''
        Gives the hull of facility j in this node as [points, not_used_points].
        :param node: The node that we are in.
        :param j: The facility.
        :param bounds: (x_lb, x_ub, y_lb, y_ub) of the facility in this node.
        :param hull_deriving_mode: 0 is that we don't update hull, 1 is cutting the root hull at all of the bounds, 2 is cutting the hull of the parent node at the last cut.
        :param h: The index of the root hull of the facility.
        '''
        root = self.root_hulls[h]
        if hull_deriving_mode == 0:
            return root

        if (node.getNumber(), j) in self.node_hulls:
            return self.node_hulls[(node.getNumber(), j)]

        hull = None
        if hull_deriving_mode == 2 and node.getParent() is not None:
            parent_hull = self.node_hulls.get((node.getParent().getNumber(), j))
            if parent_hull is not None:
                hull = parent_hull
                cut = self.get_cutval(node)
                # Only a cut over this facility changes its hull
                if cut is not None and cut[0] == j:
                    hull = self.cut_hull(parent_hull, cut[1], cut[2], cut[3])

        if hull is None:
            # We start from the root hull and cut at every bound, a bound that does not lie inside of the hull does nothing
            hull = root
            for xy in range(2):
                hull = self.cut_hull(hull, xy, bounds[2 * xy], 0)
                hull = self.cut_hull(hull, xy, bounds[2 * xy + 1], 1)

        self.node_hulls[(node.getNumber(), j)] = hull
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
            points, not_used = self.get_hull(node, j, box, self.hull_deriving_mode, h)
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

