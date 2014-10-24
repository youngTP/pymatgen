# coding: utf-8

from pymatgen.core.pbc import get_distance_and_image, cget_distance_and_image
import numpy as np

lattice = np.eye(3) * 3

def cfunc():
    cget_distance_and_image(lattice, [0.1, 0.1, 0.1], [0.9, 0.9, 0.9])

def pyfunc():
    get_distance_and_image(lattice, [0.1, 0.1, 0.1], [0.9, 0.9, 0.9])

if __name__ == "__main__":
    import timeit
    print(timeit.timeit("cfunc()", setup="from __main__ import cfunc", number=1))
    print(timeit.timeit("pyfunc()", setup="from __main__ import pyfunc",
        number=1))
