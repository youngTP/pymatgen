# coding: utf-8

from __future__ import division, unicode_literals

import itertools

import numpy as np
cimport numpy as np


DTYPE = np.float

ctypedef np.float_t DTYPE_t


cpdef cget_distance_and_image(np.ndarray[DTYPE_t, ndim=2] lattice,
                              np.ndarray[DTYPE_t, ndim=1] frac_coords1,
                              np.ndarray[DTYPE_t, ndim=1] frac_coords2,
                              jimage=None):
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
    cdef float i, j, k
    cdef np.ndarray[DTYPE_t, ndim=1] vec, coord1, coord2, image, shortest_image
    cdef float d
    cdef float dist = float('inf')
    cdef np.ndarray[DTYPE_t, ndim=1] adj1 = np.floor(frac_coords1)
    cdef np.ndarray[DTYPE_t, ndim=1] adj2 = np.floor(frac_coords2)

    if jimage is None:
        #Shift coords to unitcell
        coord1 = frac_coords1 - adj1
        coord2 = frac_coords2 - adj2
        # Generate set of images required for testing.
        # This is a cheat to create an 8x3 array of all length 3
        # combinations of 0,1
        r = (-1.0, 0.0, 1.0)
        for i, j, k in itertools.product(r, r, r):
            image = np.array([i, j, k])
            vec = coord2 - coord1 + image
            vec = np.dot(vec, lattice)
            d = (vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2) ** 0.5
            if d < dist:
                dist = d
                shortest_image = image
        return dist, adj1 - adj2 + shortest_image

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


cpdef np.ndarray[DTYPE_t, ndim=2] cget_lll_reduced_lattice(
        np.ndarray[DTYPE_t, ndim=2] lattice, float delta=0.75):
    """
    Performs a Lenstra-Lenstra-Lovasz lattice basis reduction to obtain a
    c-reduced basis. This method returns a basis which is as "good" as
    possible, with "good" defined by orthongonality of the lattice vectors.

    Args:
        delta (float): Reduction parameter. Default of 0.75 is usually
            fine.

    Returns:
        Reduced lattice.
    """
    cdef int i, k
    cdef float q

    # Transpose the lattice matrix first so that basis vectors are columns.
    # Makes life easier.
    cdef np.ndarray[DTYPE_t, ndim=2] a = lattice.T

    cdef np.ndarray[DTYPE_t, ndim=2] b = np.zeros((3, 3))  # Vectors after the
    # Gram-Schmidt process
    cdef np.ndarray[DTYPE_t, ndim=2] u = np.zeros((3, 3))  # Gram-Schmidt
    # coeffieicnts
    cdef np.ndarray[DTYPE_t, ndim=1] m = np.zeros(3)  # These are the norm
    # squared
    #  of each vec.
    cdef np.ndarray[DTYPE_t, ndim=2] p, qq, result
    cdef np.ndarray[DTYPE_t, ndim=1] v, uu

    dot = np.dot

    b[:, 0] = a[:, 0]
    m[0] = dot(b[:, 0], b[:, 0])
    for i in range(1, 3):
        u[i, 0:i] = dot(a[:, i].T, b[:, 0:i]) / m[0:i]
        b[:, i] = a[:, i] - dot(b[:, 0:i], u[i, 0:i].T)
        m[i] = dot(b[:, i], b[:, i])

    k = 2

    while k <= 3:
        # Size reduction.
        for i in range(k - 1, 0, -1):
            q = round(u[k - 1, i - 1])
            if q != 0:
                # Reduce the k-th basis vector.
                a[:, k - 1] = a[:, k - 1] - q * a[:, i - 1]
                uu = u[i - 1, 0:(i - 1)]
                uu = np.append(uu, 1)
                # Update the GS coefficients.
                u[k - 1, 0:i] = u[k - 1, 0:i] - q * uu

        # Check the Lovasz condition.
        if dot(b[:, k - 1], b[:, k - 1]) >=\
                (delta - abs(u[k - 1, k - 2]) ** 2) *\
                dot(b[:, (k - 2)], b[:, (k - 2)]):
            # Increment k if the Lovasz condition holds.
            k += 1
        else:
            #If the Lovasz condition fails,
            #swap the k-th and (k-1)-th basis vector
            v = a[:, k - 1].copy()
            a[:, k - 1] = a[:, k - 2].copy()
            a[:, k - 2] = v
            #Update the Gram-Schmidt coefficients
            for s in range(k - 1, k + 1):
                u[s - 1, 0:(s - 1)] = dot(a[:, s - 1].T,
                                          b[:, 0:(s - 1)]) / m[0:(s - 1)]
                b[:, s - 1] = a[:, s - 1] - dot(b[:, 0:(s - 1)],
                                                u[s - 1, 0:(s - 1)].T)
                m[s - 1] = dot(b[:, s - 1], b[:, s - 1])

            if k > 2:
                k -= 1
            else:
                # We have to do p/q, so do lstsq(q.T, p.T).T instead.
                p = dot(a[:, k:3].T, b[:, (k - 2):k])
                qq = np.diag(m[(k - 2):k])
                result = np.linalg.lstsq(qq.T, p.T)[0].T
                u[k:3, (k - 2):k] = result

    return a.T


def get_lll_reduced_lattice(lattice, float delta=0.75):
    """
    Performs a Lenstra-Lenstra-Lovasz lattice basis reduction to obtain a
    c-reduced basis. This method returns a basis which is as "good" as
    possible, with "good" defined by orthongonality of the lattice vectors.

    Args:
        delta (float): Reduction parameter. Default of 0.75 is usually
            fine.

    Returns:
        Reduced lattice.
    """
    # Transpose the lattice matrix first so that basis vectors are columns.
    # Makes life easier.
    a = lattice.T

    b = np.zeros((3, 3))  # Vectors after the Gram-Schmidt process
    u = np.zeros((3, 3))  # Gram-Schmidt coeffieicnts
    m = np.zeros(3)  # These are the norm squared of each vec.

    b[:, 0] = a[:, 0]


    dot = np.dot
    m[0] = dot(b[:, 0], b[:, 0])
    for i in range(1, 3):
        u[i, 0:i] = dot(a[:, i].T, b[:, 0:i]) / m[0:i]
        b[:, i] = a[:, i] - dot(b[:, 0:i], u[i, 0:i].T)
        m[i] = dot(b[:, i], b[:, i])

    k = 2

    while k <= 3:
        # Size reduction.
        for i in range(k - 1, 0, -1):
            q = round(u[k - 1, i - 1])
            if q != 0:
                # Reduce the k-th basis vector.
                a[:, k - 1] = a[:, k - 1] - q * a[:, i - 1]
                uu = list(u[i - 1, 0:(i - 1)])
                uu.append(1)
                # Update the GS coefficients.
                u[k - 1, 0:i] = u[k - 1, 0:i] - q * np.array(uu)

        # Check the Lovasz condition.
        if dot(b[:, k - 1], b[:, k - 1]) >=\
                (delta - abs(u[k - 1, k - 2]) ** 2) *\
                dot(b[:, (k - 2)], b[:, (k - 2)]):
            # Increment k if the Lovasz condition holds.
            k += 1
        else:
            #If the Lovasz condition fails,
            #swap the k-th and (k-1)-th basis vector
            v = a[:, k - 1].copy()
            a[:, k - 1] = a[:, k - 2].copy()
            a[:, k - 2] = v
            #Update the Gram-Schmidt coefficients
            for s in range(k - 1, k + 1):
                u[s - 1, 0:(s - 1)] = dot(a[:, s - 1].T,
                                          b[:, 0:(s - 1)]) / m[0:(s - 1)]
                b[:, s - 1] = a[:, s - 1] - dot(b[:, 0:(s - 1)],
                                                u[s - 1, 0:(s - 1)].T)
                m[s - 1] = dot(b[:, s - 1], b[:, s - 1])

            if k > 2:
                k -= 1
            else:
                # We have to do p/q, so do lstsq(q.T, p.T).T instead.
                p = dot(a[:, k:3].T, b[:, (k - 2):k])
                q = np.diag(m[(k - 2):k])
                result = np.linalg.lstsq(q.T, p.T)[0].T
                u[k:3, (k - 2):k] = result

    return a.T