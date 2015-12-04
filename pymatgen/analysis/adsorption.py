# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
from __future__ import absolute_import

"""
This module provides classes used to enumerate surface sites
and to find adsorption sites on slabs
"""

import numpy as np
from six.moves import range
from pymatgen.core.structure import Structure
import subprocess
import itertools
from pyhull.delaunay import DelaunayTri
from pyhull.voronoi import VoronoiTess
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

__author__ = "Joseph Montoya"
__copyright__ = "Copyright 2015, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Joseph Montoya"
__email__ = "montoyjh@lbl.gov"
__status__ = "Development"
__date__ = "December 2, 2015"


class AdsorbateSiteFinder(object):
    """
    This class finds adsorbate sites on slabs
    """

    def __init__(self, slab):
        """
        Create an AdsorbateSiteFinder object.

        Args:
            slab (Slab): slab object for which to find adsorbate
            sites
        """
        self.slab = slab

    def find_surface_sites_by_height(self, slab, window = 1.0):
        """
        This method finds surface sites by determining which sites are within
        a threshold value in height from the topmost site in a list of sites

        Args:
            site_list (list): list of sites from which to select surface sites
            window (float): threshold in angstroms of distance from topmost
                site in slab along the slab c-vector to include in surface 
                site determination
        """
        # Determine the window threshold in fractional coordinates
        c_window = window / np.linalg.norm(slab.lattice.matrix[-1])
        highest_site_z = max([site.frac_coords[-1] for site in slab.sites])
        return [site for site in slab.sites 
                if site.frac_coords[-1] >= highest_site_z - c_window]

    def find_surface_sites_by_alpha(self, slab, alpha = None):
        """
        This method finds surface sites by determining which sites are on the
        top layer of an alpha shape corresponding to the slab repeated once
        in each direction

        Args:
            site_list (list): list of sites from which to select surface sites
            alpha (float): alpha value criteria for creating alpha shape 
                for the slab object
        """
        # construct a mesh from slab repeated three times
        frac_coords = np.array([site.frac_coords for site in slab.sites])
        repeated = np.array([i + (0,) for i in 
                             itertools.product([-1,0,1], repeat=2)])
        mesh = [r + fc for r, fc in itertools.product(repeated,
                                                      frac_coords)]
        # convert mesh to input string for Clarkson hull
        mesh_string = '\n'.join([' '.join([str(j) for j in frac_coords]) 
                                 for frac_coords in mesh])
        ahull_string = subprocess.check_output(["hull", "-A"], stdin = mesh_string)
        #import pdb; pdb.set_trace()

    def get_extended_surface_mesh(self, slab, radius = 4.0, window = 1.0):
        """
        """
        surf_sites = self.find_surface_sites_by_height(slab, window = window)
        surf_str = Structure.from_sites(surf_sites)
        surface_mesh = []
        for site in surf_str.sites:
            surface_mesh += [site]
            surface_mesh += [s[0] for s in surf_str.get_neighbors(site, radius)]
            #import pdb; pdb.set_trace()
        return list(set(surface_mesh))
        #surface_mesh = list(set(surface_mesh))
        #return [site.frac_coords for site in surface_mesh]

    def find_adsorption_sites(self, slab, distance=1.0, 
                              symm_reduce=True, near_reduce=True,
                              near_reduce_threshold = 0.2):
        """
        """
        # Find vector for distance normal to x-y plane
        a, b, c = slab.lattice.matrix
        dist_vec = np.cross(a, b)
        dist_vec = distance * dist_vec / np.linalg.norm(dist_vec)
        if np.dot(dist_vec, c) < 0:
            dist_vec = -dist_vec
        # find on-top sites
        surf_sites = self.find_surface_sites_by_height(slab)
        # Get bridge sites via DelaunayTri of extended surface mesh
        mesh = self.get_extended_surface_mesh(slab)
        import pdb; pdb.set_trace()
        dt = DelaunayTri([m.frac_coords[:2] for m in mesh])

        for v in dt.vertices:
            for data in itertools.combinations(v, 2):
                bridge_sites += self.ensemble_center
        import pdb; pdb.set_trace()
        if near_reduce:
            coords = self.near_reduce(coords, threshold=near_reduce_threshold)
        if symm_reduce:
            coords = self.symm_reduce(coords)

        return coords

    def symm_reduce(self, coords_set, slab, cartesian = False):
        """

        """
        surf_sg = SpacegroupAnalyzer(slab, 0.1)
        symm_ops = surf_sg.get_symmetry_operations(cartesian = cartesian)
        full_symm_ops = generate_full_symmops(symm_ops)
        unique_coords = []
        for coords in coords_set:
            for op in full_symm_ops:
                if in_coord_list(unique_coords, op.operate(coords)):
                    unique_coords += [coords]
        return unique_pos

    def near_reduce(self, coords_set, threshold = 0.1):
        """

        """
        unique_coords = []
        for coord in coords_set:
            if in_coord_list(unique_coords, coord, threshold):
                unique_coords += [coord]
        return unique_coords

    def ensemble_center(self, site_list, indices, cartesian = False):
        """

        """
        if cartesian:
            return np.average([site_list[i].coords for i in indices], 
                              axis = 0)
        else:
            return np.average([site_list[i].frac_coords for i in indices], 
                              axis = 0)

if __name__ == "__main__":
    from pymatgen.matproj.rest import MPRester
    from pymatgen.core.surface import generate_all_slabs
    mpr = MPRester()
    struct = mpr.get_structures('Cu')[0]
    slabs = generate_all_slabs(struct, 1, 5.0, 5.0)
    asf = AdsorbateSiteFinder(slabs[0])
    #surf_sites_height = asf.find_surface_sites_by_height(slabs[0])
    #surf_sites_alpha = asf.find_surface_sites_by_alpha(slabs[0])
    sites = asf.find_adsorption_sites(slabs[0])

