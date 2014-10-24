# coding: utf-8

from pymatgen.core.pbc import get_distance_and_image, \
    cget_distance_and_image, cget_lll_reduced_lattice, get_lll_reduced_lattice
import numpy as np

lattice = np.array([[1.0, 1, 1], [-1.0, 0, 2], [3.0, 5, 6]])

def cfunc():
    # cget_distance_and_image(lattice,
    #                         np.array([0.1, 0.1, 0.1]),
    #                         np.array([0.9, 0.9, 0.9]))
    cget_lll_reduced_lattice(lattice)

def pyfunc():
    #get_distance_and_image(lattice, [0.1, 0.1, 0.1], [0.9, 0.9, 0.9])
    get_lll_reduced_lattice(lattice)

if __name__ == "__main__":
    import timeit
    print(timeit.timeit("cfunc()", setup="from __main__ import cfunc",
                        number=10))
    print(timeit.timeit("pyfunc()", setup="from __main__ import pyfunc",
        number=10))
    import cProfile, pstats, os
    cProfile.run('cfunc()', 'stats')
    p = pstats.Stats('stats')
    p.sort_stats('cumulative').print_stats(20)
    os.remove('stats')
