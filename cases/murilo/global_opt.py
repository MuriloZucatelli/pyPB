from scipy.optimize import basinhopping, shgo
from numpy import sin, sqrt, float64, infty, arange
import time
import numba
#import pyinstrument
import code
import pstats
import cProfile
import multiprocessing
import threading

a = [2, 2]


def picabol(x):
    return eggholder_jit(x, a)


#@numba.njit(parallel=True, nogil=True)
def eggholder_jit(x: float64, a):
    y = -(x[1] + 47.0) * sin(sqrt(abs(x[0] / 2.0 + (x[1] + 47.0)))) - x[0] * sin(
        sqrt(abs(x[0] - (x[1] + 47.0)))
    )

    j = infty
    for i in a:
        i + i
    for i in numba.prange(10000):
        if i < j:
            j = i
    return y


def main():

    def calbeque(xk):
        print(xk)

    def eggholder(x: float64, a):
        y = -(x[1] + 47.0) * sin(sqrt(abs(x[0] / 2.0 + (x[1] + 47.0)))) - x[0] * sin(
            sqrt(abs(x[0] - (x[1] + 47.0)))
        )
        j = infty
        for i in a:
            i + i
        for i in arange(10000):
            if i < j:
                j = i
        return y

    bounds = [(-512, 512), (-512, 512)]

    t1 = time.time()
    opt = shgo(
        lambda x: eggholder(x, a),
        bounds=bounds,
        # callback=calbeque,
        sampling_method="sobol",
        n=200,
        iters=5,
        workers=1,
    )
    t2 = time.time()
    opt_jit = shgo(
        picabol,
        bounds=bounds,
        # callback=calbeque,
        sampling_method="sobol",
        n=200,
        iters=5,
        workers=3,
    )
    t3 = time.time()

    print("normal ", opt.x, opt.fun, t2 - t1)
    print("com jit", opt_jit.x, opt_jit.fun, t3 - t2)
    print("\n", opt.nfev, opt.nlfev, opt.nit)
    print("\n", opt_jit.nfev, opt_jit.nlfev, opt_jit.nit)


if __name__ == "__main__":  # Here
    main()
