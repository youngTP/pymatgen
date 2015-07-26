__author__ = 'vivid0036'


from pymatgen.core.structure import SiteCollection, Structure
from pymatgen.core.surface import Slab
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.operations import SymmOp
from pymatgen.util.coord_utils import in_coord_list

import math
import itertools
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np
from pyhull.convex_hull import ConvexHull



"""
Wulff Shape:

1. get lattice from structure
2. get all miller index
3. combine lattice & miller index to get surfaces
4. get cross points between surfaces
"""

c = ['b', 'g', 'r', 'm', 'c', 'y']

"""

        b: blue
        g: green
        r: red
        c: cyan
        m: magenta
        y: yellow
        k: black
        w: white

"""

class wulff_3d(object):

    def __init__(self, structure, miller_list, e_surf_list):

        symmops = SpacegroupAnalyzer(structure).get_point_group_operations()

        inv_latt_matrix = structure.lattice.inv_matrix
        latt = structure.lattice
        dimension = len(miller_list[0])

        self.structure = structure
        self.unique_miller = miller_list
        self.e_surf_list = e_surf_list
        self.latt = latt
        self.inv_latt_matrix = inv_latt_matrix
        self.recp = latt.reciprocal_lattice_crystallographic
        self.symmops = symmops
        self.dimension = dimension

        normal_e_m = self.get_all_miller_e()
        # [normal, e_surf_list[i], normal_pt, dual_pt, miller]
        self.normal_e_m = normal_e_m
        print len(normal_e_m)

        normal_pt  = [x[2] for x in normal_e_m]
        dual_pt = [x[3] for x in normal_e_m]
        color_plane = [x[-2] for x in normal_e_m]


        dual_convex = ConvexHull(dual_pt)
        dual_cv_vert = dual_convex.vertices
        self.dual_cv_vert = dual_cv_vert
        # print dual_cv_vert

        # recalculate the dual of dual, get the wulff shape.
        # conner <-> surface

        wulff_pt_list = []
        wulff_pt_plane_list = []
        # dual_cv_vert
        wulff_plane_list = [[x[-1], [], [], x[-2], []] for x in normal_e_m]

        for vertices in dual_cv_vert:
            # print vertices
            i = vertices[0]
            j = vertices[1]
            k = vertices[2]
            wulff_pt = self.get_cross_pt(i, j, k)
            wulff_plane_list[i][1].append(wulff_pt)
            wulff_plane_list[j][1].append(wulff_pt)
            wulff_plane_list[k][1].append(wulff_pt)
            wulff_pt_list.append(wulff_pt)
        self.wulff_pt_list = wulff_pt_list

        wulff_convex = ConvexHull(wulff_pt_list)
        wulff_vertices = wulff_convex.vertices

        self.wulff_convex = wulff_convex
        self.wulff_vertices = wulff_vertices
        # print dual_convex.simplices


        for vertices in wulff_vertices:
            vertices.sort()
            i = vertices[0]
            j = vertices[1]
            k = vertices[2]
            v123 = dual_cv_vert[i] + dual_cv_vert[j] + dual_cv_vert[k]
            # find the plane
            for v in v123:
                if v123.count(v) == 3:
                    plane_num = v
            vertices_list = [[i, j], [j, k], [i, k]]
            wulff_plane_list[plane_num][-1] += vertices_list

            wulff_plane_list[plane_num][2].append(vertices)

        off_wulff = []
        on_wulff = []
        for plane in wulff_plane_list:
            if len(plane[-1]):
                print plane[0], 'surface on the wulff shape: color:', plane[-2]
                plane[-1].sort()
                outer_lines = []
                for line in plane[-1]:
                    if plane[-1].count(line) == 1:
                        outer_lines.append(line)
                on_wulff.append([plane[0], plane[2], plane[-2], outer_lines])
            else:
                print plane[0], 'surface not on the wulff shape'
                off_wulff.append(plane[0])
        on_wulff.sort(key= lambda x: x[-1][0], reverse=False)
        self.wulff_plane_list = wulff_plane_list
        self.on_wulff = on_wulff
        self.off_wulff = off_wulff

    def get_all_miller_e(self):
        """
        from self:
            get miller_list(unique_miller), e_surf_list and symmetry
            operations(symmops) according to lattice
        apply symmops to get all the miller index, then get normal,
        get all the planes functions for wulff shape calculation:
            |normal| = 1, e_surf is plane's distance to (0, 0, 0),
            normal[0]x + normal[1]y + normal[2]z = e_surf

        return:
            normal_e_m, item: [normal, e_surf, normal_pt, dual_pt, miller]
        """
        unique_miller = self.unique_miller
        e_surf_list = self.e_surf_list
        symmops = self.symmops
        # inv_latt_matrix = self.inv_latt_matrix
        recp = self.recp
        normal_e_m = []
        color = []
        for i in xrange(len(unique_miller)):
            color.append(c[divmod(i, 6)[1]])
        for i in xrange(len(unique_miller)):
            for op in symmops:
                miller = list(op.operate(unique_miller[i]))
                miller = [int(x) for x in miller]
                if in_coord_list(unique_miller, miller):
                    continue
                else:
                    unique_miller.append(miller)
                    e_surf_list.append(e_surf_list[i])
                    color.append(c[divmod(i, 6)[1]])
        for i in xrange(len(unique_miller)):
            miller = unique_miller[i]
            normal = recp.get_cartesian_coords(miller)
            normal /= np.linalg.norm(normal)
            e_surf = e_surf_list[i]
            normal_pt = [x*e_surf for x in normal]
            dual_pt = [x/e_surf for x in normal]
            color_plane = color[i]
            normal_e_m.append([normal, e_surf, normal_pt, dual_pt, color_plane, miller])
        # sorted by e_surf
        normal_e_m.sort(key= lambda x: x[1])

        return normal_e_m

    def get_cross_pt(self, i, j, k):
        """
        from self:
            get normal_e_m to get the plane functions
        i, j, k: plane index(same order in normal_e_m)
        """
        normal_e_m = self.normal_e_m
        matrix_surfs = [normal_e_m[i][0], normal_e_m[j][0], normal_e_m[k][0]]
        matrix_e = [normal_e_m[i][1], normal_e_m[j][1], normal_e_m[k][1]]
        cross_pt = np.dot(np.linalg.inv(matrix_surfs), matrix_e)

        return cross_pt

    def plot_wulff_pts(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for pt in self.wulff_pt_list:
            ax.scatter(pt[0], pt[1], pt[2])

        return plt
        #plt.show()


    def plot_wulff_line(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        wulff_pt_list = self.wulff_pt_list
        for plane in self.on_wulff:
            for line in plane[-1]:
                edge = [wulff_pt_list[line[0]], wulff_pt_list[line[1]]]
                ax.plot([x[0] for x in edge], [x[1] for x in edge], [x[2] for x in edge], 'k', lw=1)
        return plt
        #plt.show()

    def plot_wulff_color(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        wulff_pt_list = self.wulff_pt_list
        on_wulff = self.on_wulff

        for plane in on_wulff:
            plane_color = plane[-2]
            print plane_color, plane[0]
            pts = []
            for vertices in plane[1]:
                i = vertices[0]
                j = vertices[1]
                k = vertices[2]

                pts_vertices = [wulff_pt_list[i], wulff_pt_list[j], wulff_pt_list[k]]
                for n in range(1, 200):
                    pt1 = wulff_pt_list[i] + 0.005 * n * (wulff_pt_list[j]-wulff_pt_list[i])
                    pt2 = wulff_pt_list[i] + 0.005 * n * (wulff_pt_list[k]-wulff_pt_list[i])
                    pts_vertices.append(pt1)
                    pts_vertices.append(pt2)
                    for s in range(1, 100):
                        pt3 = (pt1*s+pt2*(100-s)) * 0.01
                        pts_vertices.append(pt3)
                pts += pts_vertices
            ax.plot([x[0] for x in pts], [x[1] for x in pts], [x[2] for x in pts], plane_color)

            for line in plane[-1]:
                edge = [wulff_pt_list[line[0]], wulff_pt_list[line[1]]]
                ax.plot([x[0] for x in edge], [x[1] for x in edge], [x[2] for x in edge], 'k', lw=1)


        plt.show()
        return plt

