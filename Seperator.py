from pyscipopt import Sepa, SCIP_RESULT

class HullSeperator(Sepa):
    def init(self, var_x, var_y, circles, A, b):
        self.var_x, self.var_y, self.circles, self.A, self.b = var_x, var_y, circles, A, b

    def sepaexeclp(self):
        model = self.model
        for i in range(len(self.A)):
            for j in range(len(self.var_x)):
                if self.A[0] * self.var_x[j] + self.A[1] * self.var_y[j] > self.b:
