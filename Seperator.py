from pyscipopt import Sepa, SCIP_RESULT
import Convex_hull as IO
import Creating_LP as CLP

def box_key(x_lb, x_ub, y_lb, y_ub, digits=9):
    return (round(x_lb, digits), round(x_ub, digits), round(y_lb, digits), round(y_ub, digits))

def facet_key(a : float,b : float,c : float, rounding : int=9):
    return (round(a,rounding), round(b,rounding), round(c,rounding))


class HullSeperator(Sepa):
    def __init__(self, var_x, var_y, circles, A, b, mode, LP_or_GS):
        self.var_x, self.var_y, self.circles, self.A, self.b, self.mode, self.LP_or_GS = var_x, var_y, circles, A, b, mode, LP_or_GS

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


    def sepaexeclp(self):
        model = self.model
        found = False 
        
        # When mode == 1 we use the seperator, otherwise the seperator gets used to update the hull after a cut
        if self.mode == 1:
            for j in range(len(self.var_x)):
                x = model.getTransformedVar(self.var_x[j])
                y = model.getTransformedVar(self.var_y[j])
                px, py = x.getLPSol(), y.getLPSol()
                # We first search if we already have a hull with the exact bounds on x and y
                key = box_key(x.getLbLocal(), x.getUbLocal(), y.getLbLocal(), y.getUbLocal())

                # We track how much each constraint is violated
                violated = []
                for i in range(len(self.A)):
                    violation = self.A[i][0] * px + self.A[i][1] * py - self.b[i]
                    if violation > 1e-6:
                        violated.append((violation, i))

                # We now pick the 5 most violated ones
                violated.sort(reverse=True)
                for violation, i in violated[:5]:
                    row = model.createEmptyRowSepa(self, name=f"hull_{j}_{i}", lhs=None, rhs=self.b[i], local=True)
                    model.cacheRowExtensions(row)
                    model.addVarToRow(row, x, self.A[i][0])
                    model.addVarToRow(row, y, self.A[i][1])
                    model.flushRowExtensions(row)
                    infeasible = model.addCut(row)
                    model.releaseRow(row)
                    if infeasible:
                        return {"result": SCIP_RESULT.CUTOFF}
                    found = True

            if found:
                return {"result": SCIP_RESULT.SEPARATED}
            return {"result": SCIP_RESULT.DIDNOTFIND}
        
        else: #self.mode == 0
            found = False
            for j in range(len(self.var_x)):        
                x = model.getTransformedVar(self.var_x[j])
                y = model.getTransformedVar(self.var_y[j])
                px, py = x.getLPSol(), y.getLPSol()

                # We add each constraint that is new
                for i in range(len(self.A)):
                    row = model.createEmptyRowSepa(self, name=f"hull_{j}_{i}", lhs=None, rhs=self.b[i], local=True)
                    model.cacheRowExtensions(row)
                    model.addVarToRow(row, x, self.A[i][0])
                    model.addVarToRow(row, y, self.A[i][1])
                    model.flushRowExtensions(row)
                    infeasible = model.addCut(row)
                    model.releaseRow(row)
                    if infeasible:
                        return {"result": SCIP_RESULT.CUTOFF}
                    found = True   
                
                if found:
                    return {"result": SCIP_RESULT.SEPARATED}
                return {"result": SCIP_RESULT.DIDNOTFIND}             

