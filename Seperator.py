from pyscipopt import Sepa, SCIP_RESULT
import Convex_hull as IO
import Creating_LP as CLP

def box_key(x_lb, x_ub, y_lb, y_ub, digits=9):
    return (round(x_lb, digits), round(x_ub, digits), round(y_lb, digits), round(y_ub, digits))

def facet_key(a : float,b : float,c : float, rounding : int=9):
    return (round(a,rounding), round(b,rounding), round(c,rounding))


class HullSeperator(Sepa):
    def __init__(self, var_x, var_y, circles, root_hulls, sep_mode=1, update=0, hull_of_facility=None):
        self.var_x, self.var_y, self.circles = var_x, var_y, circles, A, b
        self.sep_mode, self.update = sep_mode, update
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


        
    def get_hull(self, x_lb : float, x_ub : float, y_lb : float, y_ub : float, hull_deriving_mode : int, sub_hull : int, other_feasible_p : list, sorted_points : list):
        '''
        :param hull_deriving_mode: 0 is that we don't update hull, 1 is making hull from scratch, 2 is updating from earlier
        '''
        box = (x_lb, x_ub, y_lb, y_ub)
        root = self.root_hulls[sub_hull]
        if hull_deriving_mode == 0 or box == self.root_box:
            return root

        elif hull_deriving_mode == 1: # We make a new hull, using the old feasible points
            all_p = sorted_points + other_feasible_p
            feasible_p = []
            for p in all_p:
                if (x_lb <= p[0] <= x_ub) and (y_lb <= p[1] <= y_ub):
                    feasible_p.append(p)
            new_hull = IO.create_conv_region_GS(feasible_p)
            return new_hull

        else: # Mode == 2
            key = (sub_hull,) + box_key(x_lb, x_ub, y_lb, y_ub)

        return None

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
            px, py = x.getLPSol(), y.getLPSol()
            
            h = self.hull_of_facility[j] # starting hull of facility
            A, b = get_hull()

        
        # When mode == 1 we use the seperator, otherwise the seperator gets used to update the hull after a cut
        if self.sep_mode == 1:
            

                
                # We track how much each constraint is violated
                violated = []
                for i in range(len(A)):
                    violation = A[i][0] * px + A[i][1] * py - b[i]
                    if violation > 1e-6:
                        violated.append((violation, i))

                # We now pick the 5 most violated ones
                violated.sort(reverse=True)
                for violation, i in violated[:5]:
                    row = model.createEmptyRowSepa(self, name=f"hull_{j}_{i}", lhs=None, rhs=b[i], local=True)
                    model.cacheRowExtensions(row)
                    model.addVarToRow(row, x, A[i][0])
                    model.addVarToRow(row, y, A[i][1])
                    model.flushRowExtensions(row)
                    infeasible = model.addCut(row)
                    model.releaseRow(row)
                    if infeasible:
                        return {"result": SCIP_RESULT.CUTOFF}
                    found = True

            if found:
                return {"result": SCIP_RESULT.SEPARATED}
            return {"result": SCIP_RESULT.DIDNOTFIND}
        
        else: #self.sep_mode == 0
            found = False
            for j in range(len(self.var_x)):        
                x = model.getTransformedVar(self.var_x[j])
                y = model.getTransformedVar(self.var_y[j])
                px, py = x.getLPSol(), y.getLPSol()

                # We add each constraint that is new
                for i in range(len(A)):
                    key = facet_key(A[i][0], A[i][1], b[i])
                    if key in self.root_facets[h]:

        
                    row = model.createEmptyRowSepa(self, name=f"hull_{j}_{i}", lhs=None, rhs=b[i], local=True)
                    model.cacheRowExtensions(row)
                    model.addVarToRow(row, x, A[i][0])
                    model.addVarToRow(row, y, A[i][1])
                    model.flushRowExtensions(row)
                    infeasible = model.addCut(row)
                    model.releaseRow(row)
                    if infeasible:
                        return {"result": SCIP_RESULT.CUTOFF}
                    found = True   
                
                if found:
                    return {"result": SCIP_RESULT.SEPARATED}
                return {"result": SCIP_RESULT.DIDNOTFIND}             

