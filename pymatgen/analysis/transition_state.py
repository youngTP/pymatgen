# coding: utf-8

from __future__ import division, unicode_literals, print_function, absolute_import

"""
Some reimplementation of Henkelman's Transition State Analysis utilities,
which are originally in Perl. Additional features beyond those offered by
Henkelman's utilities will be added.

This allows the usage and customization in Python.
"""

__author__ = 'Shyue Ping Ong'
__copyright__ = 'Copyright 2013, The Materials Virtual Lab'
__version__ = '0.1'
__maintainer__ = 'Shyue Ping Ong'
__email__ = 'ongsp@ucsd.edu'
__date__ = '6/1/15'

import os
import glob
import math
import six
from abc import ABCMeta, abstractmethod

import numpy as np
import numpy.linalg as la
from scipy.interpolate import PiecewisePolynomial
import scipy.signal
import scipy.stats
from scipy.interpolate import interp1d

from pymatgen.util.plotting_utils import get_publication_quality_plot
from pymatgen.io.vaspio import Poscar, Outcar
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import *
from pymatgen.core.periodic_table import *
from pymatgen.io.vaspio.vasp_output import VolumetricData, Chgcar


class NEBAnalysis(object):
    """
    An NEBAnalysis class.
    """

    def __init__(self, outcars, structures, interpolation_order=3):
        """
        Initializes an NEBAnalysis from Outcar and Structure objects. Use
        the static constructors, e.g., :class:`from_dir` instead if you
        prefer to have these automatically generated from a directory of NEB
        calculations.

        Args:
            outcars ([Outcar]): List of Outcar objects. Note that these have
                to be ordered from start to end along reaction coordinates.
            structures ([Structure]): List of Structures along reaction
                coordinate. Must be same length as outcar.
            interpolation_order (int): Order of polynomial to use to
                interpolate between images. Same format as order parameter in
                scipy.interplotate.PiecewisePolynomial.
        """
        if len(outcars) != len(structures):
            raise ValueError("# of Outcars must be same as # of Structures")

        # Calculate cumulative root mean square distance between structures,
        # which serves as the reaction coordinate. Note that these are
        # calculated from the final relaxed structures as the coordinates may
        # have changed from the initial interpolation.
        r = [0]
        prev = structures[0]
        for st in structures[1:]:
            dists = np.array([s2.distance(s1) for s1, s2 in zip(prev, st)])
            r.append(np.sqrt(np.sum(dists ** 2)))
            prev = st
        r = np.cumsum(r)

        energies = []
        forces = []
        for i, o in enumerate(outcars):
            o.read_neb()
            energies.append(o.data["energy"])
            if i in [0, len(outcars) - 1]:
                forces.append(0)
            else:
                forces.append(o.data["tangent_force"])
        energies = np.array(energies)
        energies -= energies[0]
        forces = np.array(forces)
        self.r = np.array(r)
        self.energies = energies
        self.forces = forces

        # We do a piecewise interpolation between the points. Each spline (
        # cubic by default) is constrained by the boundary conditions of the
        # energies and the tangent force, i.e., the derivative of
        # the energy at each pair of points.
        self.spline = PiecewisePolynomial(
            self.r, np.array([self.energies, -self.forces]).T,
            orders=interpolation_order)

    def get_extrema(self, normalize_rxn_coordinate=True):
        """
        Returns the positions of the extrema along the MEP. Both local
        minimums and maximums are returned.

        Args:
            normalize_rxn_coordinate (bool): Whether to normalize the
                reaction coordinate to between 0 and 1. Defaults to True.

        Returns:
            (min_extrema, max_extrema), where the extrema are given as
            [(x1, y1), (x2, y2), ...].
        """
        x = np.arange(0, np.max(self.r), 0.01)
        y = self.spline(x) * 1000

        scale = 1 if not normalize_rxn_coordinate else 1 / self.r[-1]
        min_extrema = []
        max_extrema = []
        for i in range(1, len(x) - 1):
            if y[i] < y[i-1] and y[i] < y[i+1]:
                min_extrema.append((x[i] * scale, y[i]))
            elif y[i] > y[i-1] and y[i] > y[i+1]:
                max_extrema.append((x[i] * scale, y[i]))
        return min_extrema, max_extrema

    def get_plot(self, normalize_rxn_coordinate=True, label_barrier=True):
        """
        Returns the NEB plot. Uses Henkelman's approach of spline fitting
        each section of the reaction path based on tangent force and energies.

        Args:
            normalize_rxn_coordinate (bool): Whether to normalize the
                reaction coordinate to between 0 and 1. Defaults to True.
            label_barrier (bool): Whether to label the maximum barrier.

        Returns:
            matplotlib.pyplot object.
        """
        plt = get_publication_quality_plot(12, 8)
        scale = 1 if not normalize_rxn_coordinate else 1 / self.r[-1]
        x = np.arange(0, np.max(self.r), 0.01)
        y = self.spline(x) * 1000
        plt.plot(self.r * scale, self.energies * 1000, 'ro',
                 x * scale, y, 'k-', linewidth=2, markersize=10)
        plt.xlabel("Reaction coordinate")
        plt.ylabel("Energy (meV)")
        plt.ylim((np.min(y) - 10, np.max(y) * 1.02 + 20))
        if label_barrier:
            data = zip(x * scale, y)
            barrier = max(data, key=lambda d: d[1])
            plt.plot([0, barrier[0]], [barrier[1], barrier[1]], 'k--')
            plt.annotate('%.0f meV' % barrier[1],
                         xy=(barrier[0] / 2, barrier[1] * 1.02),
                         xytext=(barrier[0] / 2, barrier[1] * 1.02),
                         horizontalalignment='center')
        plt.tight_layout()
        return plt

    @classmethod
    def from_dir(cls, root_dir, relaxation_dirs=None):
        """
        Initializes a NEBAnalysis object from a directory of a NEB run.
        Note that OUTCARs must be present in all image directories. For the
        terminal OUTCARs from relaxation calculations, you can specify the
        locations using relaxation_dir. If these are not specified, the code
        will attempt to look for the OUTCARs in 00 and 0n directories,
        followed by subdirs "start", "end" or "initial", "final" in the
        root_dir. These are just some typical conventions used
        preferentially in Shyue Ping's MAVRL research group. For the
        non-terminal points, the CONTCAR is read to obtain structures. For
        terminal points, the POSCAR is used. The image directories are
        assumed to be the only directories that can be resolved to integers.
        E.g., "00", "01", "02", "03", "04", "05", "06". The minimum
        sub-directory structure that can be parsed is of the following form (
        a 5-image example is shown):

        00:
        - POSCAR
        - OUTCAR
        01, 02, 03, 04, 05:
        - CONTCAR
        - OUTCAR
        06:
        - POSCAR
        - OUTCAR

        Args:
            root_dir (str): Path to the root directory of the NEB calculation.
            relaxation_dirs (tuple): This specifies the starting and ending
                relaxation directories from which the OUTCARs are read for the
                terminal points for the energies.

        Returns:
            NEBAnalysis object.
        """
        neb_dirs = []

        for d in os.listdir(root_dir):
            pth = os.path.join(root_dir, d)
            if os.path.isdir(pth) and d.isdigit():
                i = int(d)
                neb_dirs.append((i, pth))
        neb_dirs = sorted(neb_dirs, key=lambda d: d[0])
        outcars = []
        structures = []

        # Setup the search sequence for the OUTCARs for the terminal
        # directories.
        terminal_dirs = []
        if relaxation_dirs is not None:
            terminal_dirs.append(relaxation_dirs)
        terminal_dirs.append((neb_dirs[0][1], neb_dirs[-1][1]))
        terminal_dirs.append([os.path.join(root_dir, d)
                              for d in ["start", "end"]])
        terminal_dirs.append([os.path.join(root_dir, d)
                              for d in ["initial", "final"]])

        for i, d in neb_dirs:
            outcar = glob.glob(os.path.join(d, "OUTCAR*"))
            contcar = glob.glob(os.path.join(d, "CONTCAR*"))
            poscar = glob.glob(os.path.join(d, "POSCAR*"))
            terminal = i == 0 or i == neb_dirs[-1][0]
            if terminal:
                found = False
                for ds in terminal_dirs:
                    od = ds[0] if i == 0 else ds[1]
                    outcar = glob.glob(os.path.join(od, "OUTCAR*"))
                    if outcar:
                        outcar = sorted(outcar)
                        outcars.append(Outcar(outcar[-1]))
                        found = True
                        break
                if not found:
                    raise ValueError("OUTCAR cannot be found for terminal "
                                     "point %s" % d)
                structures.append(Poscar.from_file(poscar[0]).structure)
            else:
                outcars.append(Outcar(outcar[0]))
                structures.append(Poscar.from_file(contcar[0]).structure)
        return NEBAnalysis(outcars, structures)
        

class NEBPathfinder:
    def __init__(self, start_struct, end_struct, relax_sites, v, n_images=20):
        """
        General pathfinder for interpolating between two structures, where the interpolating path is calculated with
        the elastic band method with respect to the given static potential for sites whose indices are given in
        relax_sites, and is linear otherwise.
        :param start_struct, end_struct - Endpoint structures to interpolate between
        :param relax_sites - List of site indices whose interpolation paths should be relaxed
        :param v - Static potential field to use for the elastic band relaxation
        :param n_images - Number of interpolation images to generate
        """
        self.__s1 = start_struct
        self.__s2 = end_struct
        self.__relax_sites = relax_sites
        self.__v = v
        self.__n_images = n_images
        self.__images = None
        self.interpolate()

    def interpolate(self):
        """
        Finds a set of n_images from self.s1 to self.s2, where all sites except for the ones given in relax_sites,
        the interpolation is linear (as in pymatgen.core.structure.interpolate), and for the site indices given
        in relax_sites, the path is relaxed by the elastic band method within the static potential V.
        :param relax_sites: List of site indices for which the interpolation path needs to be relaxed
        :param v: static potential to use for relaxing the interpolation path
        :param n_images: number of images to generate along the interpolation path
        """
        images = self.__s1.interpolate(self.__s2, nimages=self.__n_images, interpolate_lattices=True)
        for site_i in self.__relax_sites:
            start_f = images[0].sites[site_i].frac_coords
            end_f = images[-1].sites[site_i].frac_coords

            path = NEBPathfinder.string_relax(NEBPathfinder.__f2d(start_f, self.__v),
                                              NEBPathfinder.__f2d(end_f, self.__v),
                                              self.__v, n_images=(self.__n_images+1),
                                                                    dr=[self.__s1.lattice.a/self.__v.shape[0],
                                                                        self.__s1.lattice.b/self.__v.shape[1],
                                                                        self.__s1.lattice.c/self.__v.shape[2]])
            for image_i, image in enumerate(images):
                image.translate_sites(site_i,
                                      NEBPathfinder.__d2f(path[image_i], self.__v) - image.sites[site_i].frac_coords,
                                      frac_coords=True, to_unit_cell=True)
        self.__images = images

    @property
    def images(self):
        """
        Returns a list of structures interpolating between the start and endpoint structures.
        """
        return self.__images

    def plot_images(self, outfile):
        """
        Generates a POSCAR with the calculated diffusion path with respect to the first endpoint.
        :param outfile: Output file for the POSCAR
        """
        sum_struct = self.__images[0].sites
        for image in self.__images:
            for site_i in self.__relax_sites:
                sum_struct.append(PeriodicSite(image.sites[site_i].specie, image.sites[site_i].frac_coords,
                                               self.__images[0].lattice, to_unit_cell=True, coords_are_cartesian=False))
        sum_struct = Structure.from_sites(sum_struct, validate_proximity=False)
        p = Poscar(sum_struct)
        p.write_file(outfile)

    @staticmethod
    def string_relax(start, end, V, n_images=25, dr=None, h=3.0, k=0.17, min_iter=100, max_iter=10000, max_tol=5e-6):
        """
        Implements path relaxation via the elastic band method. In general, the method is to define a path by a set of
        points (images) connected with bands with some elasticity constant k. The images then relax along the forces
        found in the potential field V, counterbalanced by the elastic response of the elastic band. In general the
        endpoints of the band can be allowed to relax also to their local minima, but in this calculation they are kept
        fixed.
        :param start, end - Endpoints of the path calculation given in discrete coordinates with respect to the grid in V
        :param V - potential field through which to calculate the path
        :param n_images - number of images used to define the path. In general anywhere from 20 to 40 seems to be good.
        :param dr - Conversion ratio from discrete coordinates to real coordinates for each of the three coordinate vectors
        :param h - Step size for the relaxation. h = 0.1 works reliably, but is slow. h=10 diverges with large gradients
                    but for the types of gradients seen in CHGCARs, works pretty reliably
        :param k - Elastic constant for the band (in real units, not discrete)
        :param min_iter, max_iter - Number of optimization steps the string will take before exiting (even if unconverged)
        :param max_tol - Convergence threshold such that if the string moves by less than max_tol in a step, and at least
                        min_iter steps have passed, the algorithm will terminate. Depends strongly on the size of the
                        gradients in V, but 5e-6 works reasonably well for CHGCARs

        """
        #
        # This code is based on the MATLAB example provided by
        # Prof. Eric Vanden-Eijnden of NYU
        # (http://www.cims.nyu.edu/~eve2/main.htm)
        #

        print("Getting path from {} to {} (coords wrt V grid)".format(start, end))
        # Set parameters
        if not dr:
            dr = np.array([1.0/V.shape[0], 1.0/V.shape[1], 1.0/V.shape[2]])
        else:
            dr = np.array(dr, dtype=float)
        keff = k * dr * n_images
        h0 = h

        # Initialize string
        g1 = np.linspace(0, 1, n_images)
        s0 = start
        s1 = end
        s = np.array([g * (s1-s0) for g in g1]) + s0
        ds = s - np.roll(s,1,axis=0)
        ds[0] = (ds[0] - ds[0])
        ls = np.cumsum(la.norm(ds, axis=1))
        ls = ls/ls[-1]
        fi = interp1d(ls,s,axis=0)
        s = fi(g1)

        # Evaluate initial distances (for elastic equilibrium)
        ds0_plus = s - np.roll(s,1,axis=0)
        ds0_minus =  s - np.roll(s,-1,axis=0)
        ds0_plus[0] = (ds0_plus[0] - ds0_plus[0])
        ds0_minus[-1] = (ds0_minus[-1] - ds0_minus[-1])

        # Evolve string
        for step in range(0, max_iter):
            if step > min_iter:
                h = h0 * np.exp(-2.0 * (step - min_iter)/max_iter) # Gradually decay step size to prevent oscillations
            else:
                h = h0
            # Calculate forces acting on string
            dV = np.gradient(V)
            d = V.shape
            s0 = s
            edV = np.array([[dV[0][int(pt[0])%d[0]][int(pt[1])%d[1]][int(pt[2])%d[2]] / dr[0],
                             dV[1][int(pt[0])%d[0]][int(pt[1])%d[1]][int(pt[2])%d[2]] / dr[0],
                             dV[2][int(pt[0])%d[0]][int(pt[1])%d[1]][int(pt[2])%d[2]] / dr[0]] for pt in s])
            #if(step % 100 == 0):
            #    print(edV)

            # Update according to force due to potential and string elasticity
            ds_plus = s - np.roll(s,1,axis=0)
            ds_minus =  s - np.roll(s,-1,axis=0)
            ds_plus[0] = (ds_plus[0] - ds_plus[0])
            ds_minus[-1] = (ds_minus[-1] - ds_minus[-1])
            Fpot = edV
            Fel = keff * (la.norm(ds_plus) - la.norm(ds0_plus)) * (ds_plus / la.norm(ds_plus))
            Fel += keff * (la.norm(ds_minus) - la.norm(ds0_minus)) * (ds_minus / la.norm(ds_minus))
            s = s - h * (Fpot + Fel)

            # Fix endpoints
            s[0] = s0[0]
            s[-1] = s0[-1]

            # Reparametrize string
            ds = s - np.roll(s,1,axis=0)
            ds[0] = (ds[0] - ds[0])
            ls = np.cumsum(la.norm(ds, axis=1))
            ls = ls/ls[-1]
            fi = interp1d(ls,s,axis=0)
            s = fi(g1)

            tol = la.norm((s-s0) * dr) / n_images / h

            if (tol > 1e10):
                raise ValueError("Pathfinding failed, path diverged! Consider reducing h to avoid divergence.")

            if (step > min_iter and tol < max_tol):
                print("Converged at step {}".format(step))
                break

            if (step % 100 == 0):
                print ("Step {} - ds = {}".format(step, tol))
        return s

    @staticmethod
    def __f2d(frac_coords, v):
        """
        Converts fractional coordinates to discrete coordinates with respect to the grid size of v
        """
        #frac_coords = frac_coords % 1
        return np.array([int(frac_coords[0]*v.shape[0]),
                         int(frac_coords[1]*v.shape[1]),
                         int(frac_coords[2]*v.shape[2])])

    @staticmethod
    def __d2f(disc_coords, v):
        """
        Converts a point given in discrete coordinates withe respect to the grid in v to fractional coordinates.
        """
        return np.array([disc_coords[0]/v.shape[0],
                         disc_coords[1]/v.shape[1],
                         disc_coords[2]/v.shape[2]])

class StaticPotential(six.with_metaclass(ABCMeta)):
    """
    Defines a general static potential for diffusion calculations. Implements grid-rescaling and smearing for the
    potential grid. Also provides a function to normalize the potential from 0 to 1 (recommended).
    """
    def __init__(self, struct, pot):
        self.__v = pot
        self.__s = struct

    def get_v(self):
        """
        Returns the potential
        """
        return self.__v

    def normalize(self):
        """
        Sets the potential range 0 to 1.
        """
        self.__v = self.__v - np.amin(self.__v)
        self.__v = self.__v / np.amax(self.__v)

    def rescale_field(self, new_dim):
        """
        Changes the discretization of the potential field by linear interpolation. This is necessary if the potential field
        obtained from DFT is strangely skewed, or is too fine or coarse. Obeys periodic boundary conditions at the edges of
        the cell. Alternatively useful for mixing potentials that originally are on different grids.

        :param new_dim: tuple giving the numpy shape of the new grid
        """
        v_dim = self.__v.shape
        padded_v = np.lib.pad(self.__v, ((0,1), (0,1), (0,1)), mode='wrap')
        ogrid_list = np.array([list(c) for c in list(np.ndindex(v_dim[0]+1, v_dim[1]+1, v_dim[2]+1))])
        v_ogrid = padded_v.reshape(((v_dim[0]+1) * (v_dim[1]+1) * (v_dim[2]+1), -1))
        ngrid_a, ngrid_b, ngrid_c = np.mgrid[0 : v_dim[0] : v_dim[0]/new_dim[0],
                                             0 : v_dim[1] : v_dim[1]/new_dim[1],
                                             0 : v_dim[2] : v_dim[2]/new_dim[2]]

        v_ngrid = scipy.interpolate.griddata(ogrid_list, v_ogrid, (ngrid_a, ngrid_b, ngrid_c), method='linear').reshape((new_dim[0], new_dim[1], new_dim[2]))
        self.__v = v_ngrid

    def gaussian_smear(self, r):
        """
        Applies an isotropic Gaussian smear of width (standard deviation) r to the potential field. This is necessary to
        avoid finding paths through narrow minima or nodes that may exist in the field (although any potential or
        charge distribution generated from GGA should be relatively smooth anyway). The smearing obeys periodic
        boundary conditions at the edges of the cell.

        :param r - Smearing width in cartesian coordinates, in the same units as the structure lattice vectors
        """
        # Since scaling factor in fractional coords is not isotropic, have to have different radii in 3 directions
        a_lat = self.__s.lattice.a
        b_lat = self.__s.lattice.b
        c_lat = self.__s.lattice.c

        # Conversion factors for discretization of v
        v_dim = self.__v.shape
        r_frac = (r / a_lat, r / b_lat, r / c_lat)
        r_disc = (int(math.ceil(r_frac[0] * v_dim[0])), int(math.ceil(r_frac[1] * v_dim[1])),
                  int(math.ceil(r_frac[2] * v_dim[2])))

        # Apply smearing
        # Gaussian filter
        gauss_dist = np.zeros((r_disc[0] * 4 + 1, r_disc[1] * 4 + 1, r_disc[2] * 4 + 1))
        for g_a in np.arange(-2.0 * r_disc[0], 2.0 * r_disc[0] + 1, 1.0):
            for g_b in np.arange(-2.0 * r_disc[1], 2.0 * r_disc[1] + 1, 1.0):
                for g_c in np.arange(-2.0 * r_disc[2], 2.0 * r_disc[2] + 1, 1.0):
                    g = np.array([g_a / v_dim[0], g_b / v_dim[1], g_c / v_dim[2]]).T
                    gauss_dist[int(g_a + r_disc[0])][int(g_b + r_disc[1])][int(g_c + r_disc[2])] = la.norm(np.dot(self.__s.lattice.matrix, g))/r
        gauss = scipy.stats.norm.pdf(gauss_dist)
        gauss = gauss/np.sum(gauss, dtype=float)
        padded_v = np.pad(self.__v, ((r_disc[0], r_disc[0]), (r_disc[1], r_disc[1]), (r_disc[2], r_disc[2])), mode='wrap')
        smeared_v = scipy.signal.convolve(padded_v, gauss, mode='valid')
        self.__v = smeared_v

class ChgcarPotential(StaticPotential):
    '''
    Implements a potential field based on the charge density output from VASP.
    '''
    def __init__(self, chgcar, smear=False, normalize=True):
        """
        :param chgcar: Chgcar object based on a VASP run of the structure of interest (Chgcar.from_file("CHGCAR"))
        :param smear: Whether or not to apply a Gaussian smearing to the potential
        :param normalize: Whether or not to normalize the potential to range from 0 to 1
        """
        v = chgcar.data['total']
        v = v / (v.shape[0] * v.shape[1] * v.shape[2])
        StaticPotential.__init__(self, chgcar.structure, v)
        if smear:
            self.gaussian_smear(2.0)
        if normalize:
            self.normalize()

class FreeVolumePotential(StaticPotential):
    '''
    Implements a potential field based on geometric distances from atoms in the structure - basically, the potential
    is lower at points farther away from any atoms in the structure.
    '''
    def __init__(self, struct, dim, smear=False, normalize=True):
        """
        :param struct: Unit cell on which to base the potential
        :param dim: Grid size for the potential
        :param smear: Whether or not to apply a Gaussian smearing to the potential
        :param normalize: Whether or not to normalize the potential to range from 0 to 1
        """
        self.__s = struct
        v = FreeVolumePotential.__add_gaussians(struct, dim)
        StaticPotential.__init__(self, struct, v)
        if smear:
            self.gaussian_smear(2.0)
        if normalize:
            self.normalize()

    @staticmethod
    def __add_gaussians(s, dim, r=1.5):
        gauss_dist = np.zeros(dim)
        for a_d in np.arange(0.0, dim[0], 1.0):
            for b_d in np.arange(0.0, dim[1], 1.0):
                for c_d in np.arange(0.0, dim[2], 1.0):
                    coords_f = np.array([a_d / dim[0], b_d / dim[1], c_d / dim[2]])
                    d_f = sorted(s.get_sites_in_sphere(coords_f, s.lattice.a), key=lambda x:x[1])[0][1]
                    #print(d_f)
                    gauss_dist[int(a_d)][int(b_d)][int(c_d)] = d_f / r
        v = scipy.stats.norm.pdf(gauss_dist)
        return v

class MixedPotential(StaticPotential):
    '''
    Implements a potential that is a weighted sum of some other potentials
    '''
    def __init__(self, potentials, coefficients, smear=False, normalize=True):
        """
        :param potentials: List of objects extending the StaticPotential superclass
        :param coefficients: Mixing weights for the elements of the potentials list
        :param smear: Whether or not to apply a Gaussian smearing to the potential
        :param normalize: Whether or not to normalize the potential to range from 0 to 1
        """
        v = potentials[0].get_v() * coefficients[0]
        s = potentials[0].__s
        for i in range(1, len(potentials)):
            v += potentials[i].get_v() * coefficients[i]
        StaticPotential.__init__(self, s, v)
        if smear:
            self.gaussian_smear(2.0)
        if normalize:
            self.normalize()
