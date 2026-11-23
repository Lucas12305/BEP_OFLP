from pyscipopt import Sepa, SCIP_RESULT

def box_key(x_lb, x_ub, y_lb, y_ub, digits=9):
    return (round(x_lb, digits), round(x_ub, digits), round(y_lb, digits), round(y_ub, digits))


class HullSeperator(Sepa):
    def __init__(self, var_x, var_y, circles, A, b, mode):
        self.var_x, self.var_y, self.circles, self.A, self.b, self.mode = var_x, var_y, circles, A, b, mode



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

