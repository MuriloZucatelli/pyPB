import numpy as np
from pymoo.core.problem import Problem, ElementwiseProblem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.termination import get_termination
from pymoo.optimize import minimize
import matplotlib.pyplot as plt


class SphereWithConstraint(Problem):

    def __init__(self):
        super().__init__(n_var=10, n_obj=1, n_ieq_constr=1, xl=0.0, xu=1.0)

    def _evaluate(self, x, out, *args, **kwargs):
        out["F"] = np.sum((x - 0.5) ** 2, axis=1)
        out["G"] = 0.1 - out["F"]


class MyProblem2(Problem):

    def __init__(self):
        super().__init__(n_var=2, n_obj=3, xl=-2.0, xu=2.0)

    def _evaluate(self, x, out, *args, **kwargs):
        f1 = 100 * (x[:, 0] ** 2 + x[:, 1] ** 2)
        f2 = (x[:, 0] - 1) ** 2 + x[:, 1] ** 2
        f3 = (x[:, 0] - 1) ** 3 + x[:, 1] ** 3
        out["F"] = np.column_stack([f1, f2, f3])


class MyProblem(ElementwiseProblem):

    def __init__(self):
        super().__init__(
            n_var=2,
            n_obj=11,
            n_ieq_constr=2,
            xl=np.array([-2, -2]),
            xu=np.array([2, 2]),
        )

    def _evaluate(self, x, out, *args, **kwargs):
        f1 = 100 * (x[0] ** 2 + x[1] ** 2)
        f2 = (x[0] - 1) ** 2 + x[1] ** 2
        f3 = (x[0] - 1) + x[1]
        f4 = (x[0] - 2) + x[1]
        f5 = (x[0] - 5) + x[1]
        f6 = (x[0] + 11) + x[1]
        f7 = (x[0] * 41) + x[1]
        f8 = (x[0] + 31) + x[1]
        f9 = (x[0] - 21) + x[1]
        f10 = (x[0] - 51) + x[1]
        f11 = (x[0] - 71) + x[1]

        g1 = 2 * (x[0] - 0.1) * (x[0] - 0.9) / 0.18
        g2 = -20 * (x[0] - 0.4) * (x[0] - 0.6) / 4.8

        out["F"] = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11]
        out["G"] = [g1, g2]


problem = MyProblem2()

algorithm = NSGA2(
    pop_size=40,
    n_offsprings=10,
    sampling=FloatRandomSampling(),
    crossover=SBX(prob=0.9, eta=15),
    mutation=PM(eta=20),
    eliminate_duplicates=True,
)

termination = get_termination("n_gen", 40)

res = minimize(problem, algorithm, termination, seed=1, save_history=True, verbose=True)

X = res.X
F = res.F

plt.figure(figsize=(7, 5))
plt.scatter(F[:, 0], F[:, 1], s=30, facecolors="none", edgecolors="blue")
plt.title("Objective Space")
plt.show()
