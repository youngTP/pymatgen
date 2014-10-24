# coding: utf-8

from __future__ import division, unicode_literals

import itertools

import numpy as np
cimport numpy as np


DTYPE = np.float

ctypedef np.float_t DTYPE_t


cpdef cget_distance_and_image(np.ndarray[DTYPE_t, ndim=2] lattice,
                             frac_coords1, frac_coords2, jimage=None):
    """
    Gets distance between two frac_coords assuming periodic boundary
    conditions. If the index jimage is not specified it selects the j
    image nearest to the i atom and returns the distance and jimage
    indices in terms of lattice vector translations. If the index jimage
    is specified it returns the distance between the frac_coords1 and
    the specified jimage of frac_coords2, and the given jimage is also
    returned.

    Args:
        fcoords1 (3x1 array): Reference fcoords to get distance from.
        fcoords2 (3x1 array): fcoords to get distance from.
        jimage (3x1 array): Specific periodic image in terms of
            lattice translations, e.g., [1,0,0] implies to take periodic
            image that is one a-lattice vector away. If jimage == None,
            the image that is nearest to the site is found.

    Returns:
        (distance, jimage): distance and periodic lattice translations
        of the other site for which the distance applies. This means that
        the distance between frac_coords1 and (jimage + frac_coords2) is
        equal to distance.
    """
    cdef np.ndarray[DTYPE_t, ndim=1] adj1, adj2, coord1, coord2, vec
    cdef float d
    cdef float dist = float('inf')

    if jimage is None:
        #The following code is heavily vectorized to maximize speed.
        #Get the image adjustment necessary to bring coords to unit_cell.
        adj1 = np.floor(frac_coords1)
        adj2 = np.floor(frac_coords2)
        #Shift coords to unitcell
        coord1 = frac_coords1 - adj1
        coord2 = frac_coords2 - adj2
        # Generate set of images required for testing.
        # This is a cheat to create an 8x3 array of all length 3
        # combinations of 0,1
        r = (-1, 0, 1)

        for i, j, k in itertools.product(r, r, r):
            vec = coord2 - coord1 + [i, j, j]
            vec = np.dot(vec, lattice)
            d = np.linalg.norm(vec)
            if d < dist:
                dist = d
                image = [i, j, k]
        return dist, adj1 - adj2 + image

    mapped_vec = np.dot(jimage + frac_coords2 - frac_coords1, lattice)
    return np.linalg.norm(mapped_vec), jimage


#Non-optimized version
def get_distance_and_image(lattice, frac_coords1, frac_coords2, jimage=None):
    """
    Gets distance between two frac_coords assuming periodic boundary
    conditions. If the index jimage is not specified it selects the j
    image nearest to the i atom and returns the distance and jimage
    indices in terms of lattice vector translations. If the index jimage
    is specified it returns the distance between the frac_coords1 and
    the specified jimage of frac_coords2, and the given jimage is also
    returned.

    Args:
        fcoords1 (3x1 array): Reference fcoords to get distance from.
        fcoords2 (3x1 array): fcoords to get distance from.
        jimage (3x1 array): Specific periodic image in terms of
            lattice translations, e.g., [1,0,0] implies to take periodic
            image that is one a-lattice vector away. If jimage == None,
            the image that is nearest to the site is found.

    Returns:
        (distance, jimage): distance and periodic lattice translations
        of the other site for which the distance applies. This means that
        the distance between frac_coords1 and (jimage + frac_coords2) is
        equal to distance.
    """
    if jimage is None:
        #The following code is heavily vectorized to maximize speed.
        #Get the image adjustment necessary to bring coords to unit_cell.
        adj1 = np.floor(frac_coords1)
        adj2 = np.floor(frac_coords2)
        #Shift coords to unitcell
        coord1 = frac_coords1 - adj1
        coord2 = frac_coords2 - adj2
        # Generate set of images required for testing.
        # This is a cheat to create an 8x3 array of all length 3
        # combinations of 0,1
        test_set = np.unpackbits(np.array([5, 57, 119],
                                          dtype=np.uint8)).reshape(8, 3)
        images = np.copysign(test_set, coord1 - coord2)
        # Create tiled cartesian coords for computing distances.
        vec = np.tile(coord2 - coord1, (8, 1)) + images
        vec = np.dot(vec, lattice)
        # Compute distances manually.
        dist = np.sqrt(np.sum(vec ** 2, 1)).tolist()
        # Return the minimum distance and the adjusted image corresponding
        # to the min distance.
        mindist = min(dist)
        ind = dist.index(mindist)
        return mindist, adj1 - adj2 + images[ind]

    mapped_vec = np.dot(jimage + frac_coords2 - frac_coords1, lattice)
    return np.linalg.norm(mapped_vec), jimage