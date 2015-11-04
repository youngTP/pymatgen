__author__ = 'vivid0036'


from pymatgen.core.structure import SiteCollection, Structure
from pymatgen.core.surface import Slab
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.operations import SymmOp
from pymatgen.util.coord_utils import in_coord_list


import math
import itertools
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
# from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cm
import numpy as np
from scipy import array
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

"""
Notes for brewer color plotting methods:
get brewer_color from (distinct values)
    http://colorbrewer2.org/
make sure the input e_surf are from small to large,
and the brewer color from dark to light
(since small e_surf tends to have larger area on wulff shapes)
    some color set of 5:
        1. Leaves: ['#bd0026', '#f03b20', '#fd8d3c', '#fecc5c', '#ffffb2']
        2. Milk-Wine: ['#b30000', '#e34a33', '#fc8d59', '#fdcc8a', '#fef0d9']
        3. Cloud-Night: ['#810f7c', '#8856a7', '#8c96c6', '#b3cde3', '#edf8fb']
        4. Cherry-Violet: ['#7a0177', '#c51b8a', '#f768a1', '#fbb4b9', '#feebe2']
        5. Boat-Ocean: ['#253494', '#2c7fb8', '#41b6c4', '#a1dab4','#a1dab4', '#ffffcc']
    But when set values to brewer_color, make sure every surface have a related color
    (len(e_surf_list)=len(brewer_color))
"""

class wulff_3d(object):

    def __init__(self, structure, miller_list, e_surf_list, bar_range=[], color_set='autumn', brewer_color=[],
                 color_shift=0, grid_off=True, axis_off=True, show_area=True, label_miller=True, alpha=0.5):

        symmops = SpacegroupAnalyzer(structure).get_point_group_operations()

        inv_latt_matrix = structure.lattice.inv_matrix
        latt = structure.lattice
        dimension = len(miller_list[0])

        # Set color mapping according to E_surf
        for x in e_surf_list:
            bar_range.append(x)
        print bar_range
        c_map = plt.get_cmap(color_set)
        cnorm = colors.Normalize(vmin=min(bar_range), vmax=max(bar_range) * (1 + color_shift * len(e_surf_list)))
        self.cnorm = cnorm
        scalar_map = cm.ScalarMappable(norm=cnorm, cmap=c_map)

        e_surf_list_cplot = []
        for i in xrange(len(e_surf_list)):
            print e_surf_list[i] + i * color_shift * (max(bar_range) - min(bar_range))
            e_surf_list_cplot.append(e_surf_list[i] + i * color_shift * (max(bar_range) - min(bar_range)))

        color_e_surf = [scalar_map.to_rgba(x) for x in e_surf_list_cplot]
        print  '^ ^', color_e_surf
        scalar_map.set_array(bar_range)

        color_proxy = [plt.Rectangle((2, 2), 1, 1, fc=x, alpha=alpha) for x in color_e_surf]
        bar_range.sort()
        self.brewer_color = brewer_color
        self.scalarcm = scalar_map
        self.color_e_surf = color_e_surf
        self.bar_range = bar_range
        self.color_proxy = color_proxy
        self.color_set = color_set
        self.label_miller = label_miller
        self.show_area = show_area
        self.alpha = alpha
        self.axis_off = axis_off
        self.grid_off = grid_off

        self.structure = structure
        self.input_miller = [str(x) for x in miller_list]
        self.input_esurf = [x for x in e_surf_list]
        self.unique_miller = miller_list
        self.e_surf_list = e_surf_list
        self.latt = latt
        self.inv_latt_matrix = inv_latt_matrix
        self.recp = latt.reciprocal_lattice_crystallographic
        self.symmops = symmops
        self.dimension = dimension

        normal_e_m = self.get_all_miller_e()
        # [normal, e_surf, normal_pt, dual_pt, color_plane, m_ind_orig, miller]
        self.normal_e_m = normal_e_m
        print len(normal_e_m)

        normal_pt  = [x[2] for x in normal_e_m]
        dual_pt = [x[3] for x in normal_e_m]
        color_plane = [x[4] for x in normal_e_m]

        dual_convex = ConvexHull(dual_pt)
        dual_cv_vert = dual_convex.vertices
        self.dual_cv_vert = dual_cv_vert
        # print dual_cv_vert

        # recalculate the dual of dual, get the wulff shape.
        # conner <-> surface

        wulff_pt_list = []
        wulff_pt_plane_list = []
        # dual_cv_vert
        # [miller, wulff_pt]
        wulff_plane_list = [[x[-1], [], [], x[4], x[5], []] for x in normal_e_m]

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
                # print plane[0], 'surface on the wulff shape: color:', plane[-2]
                plane[-1].sort()
                outer_lines = []
                for line in plane[-1]:
                    if plane[-1].count(line) == 1:
                        outer_lines.append(line)
                on_wulff.append([plane[0], plane[1], plane[2], plane[3], plane[4], outer_lines])
            else:
                # print plane[0], 'surface not on the wulff shape'
                off_wulff.append(plane[0])
        on_wulff.sort(key= lambda x: x[-1][0], reverse=False)
        self.wulff_plane_list = wulff_plane_list
        self.on_wulff = on_wulff
        self.off_wulff = off_wulff

        self.color_area_list = self.get_wulff_area()
        miller_area = []
        for m in xrange(len(self.input_miller)):
            print m
            miller_area.append(self.input_miller[m] + ' Total Areas: ' + str(round(self.color_area_list[m], 4)))
        self.miller_area = miller_area
        print miller_area
        print len(e_surf_list)

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
            normal_e_m, item: [normal, e_surf, normal_pt, dual_pt,
            color_plane, m_ind_orig, miller]
        """
        unique_miller = self.unique_miller
        unique_miller_ind = list(enumerate(unique_miller))
        e_surf_list = self.e_surf_list
        symmops = self.symmops
        # inv_latt_matrix = self.inv_latt_matrix
        recp = self.recp
        color_e_surf = self.color_e_surf
        normal_e_m = []
        color = []
        miller_ind_orig = [x[0] for x in unique_miller_ind]
        for i in xrange(len(unique_miller_ind)):
            color.append(color_e_surf[divmod(i, len(color_e_surf))[1]])
        for i in xrange(len(unique_miller_ind)):
            for op in symmops:
                miller = list(op.operate(unique_miller[i]))
                miller = [int(x) for x in miller]
                if in_coord_list(unique_miller, miller):
                    continue
                else:
                    unique_miller.append(miller)
                    e_surf_list.append(e_surf_list[i])
                    miller_ind_orig.append(i)
                    color.append(color_e_surf[divmod(i, len(color_e_surf))[1]])

        for i in xrange(len(unique_miller)):
            miller = unique_miller[i]
            normal = recp.get_cartesian_coords(miller)
            normal /= np.linalg.norm(normal)
            e_surf = e_surf_list[i]
            normal_pt = [x*e_surf for x in normal]
            dual_pt = [x/e_surf for x in normal]
            color_plane = color[i]
            m_ind_orig = miller_ind_orig[i]
            normal_e_m.append([normal, e_surf, normal_pt, dual_pt,
                               color_plane, m_ind_orig, miller])
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

    def get_wulff_area(self):
        wulff_pt_list = self.wulff_pt_list
        on_wulff = self.on_wulff
        color_e_surf = self.color_e_surf
        # initialize
        ind_area = [0] * len(color_e_surf)
        plane_center_label = []

        for plane in on_wulff:
            plane_color = plane[3]
            plane_m_ind_orig = plane[4]
            plane_d = plane[2]
            print plane_color, plane[0]
            plane_area = 0

            for vertices in plane[2]:
                i = vertices[0]
                j = vertices[1]
                k = vertices[2]
                pts_vertices = [wulff_pt_list[i], wulff_pt_list[j], wulff_pt_list[k]]
                v1 = pts_vertices[1] - pts_vertices[0]
                v2 = pts_vertices[2] - pts_vertices[0]
                # print v1, v2
                area_tri = np.linalg.norm(np.cross(v1, v2)) / 2
                print area_tri
                plane_area += area_tri
            print 'plane_area', plane_area

            for i in xrange(len(self.input_miller)):
                if plane_m_ind_orig == i:
                    ind_area[i] += plane_area
        print color_e_surf, ind_area

        return ind_area



    def plot_wulff_pts(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for pt in self.wulff_pt_list:
            ax.scatter(pt[0], pt[1], pt[2])
        plt.gca().set_aspect('equal', adjustable='box')
        if self.show_area == True:
            ax.legend(self.color_proxy, self.miller_area, loc='upper left',
                      bbox_to_anchor=(-0.2, 1.05), fancybox=True, shadow=True)
        else:
            ax.legend(self.color_proxy, self.input_miller, loc='upper center',
                      bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True, shadow=True)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        fig.colorbar(self.scalarcm)

        if self.grid_off == True:
            ax.grid('off')
        if self.axis_off == True:
            ax.axis('off')

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
        plt.gca().set_aspect('equal', adjustable='box')
        if self.show_area == True:
            ax.legend(self.color_proxy, self.miller_area, loc='upper left',
                      bbox_to_anchor=(-0.2, 1.05), fancybox=True, shadow=True)
        else:
            ax.legend(self.color_proxy, self.input_miller, loc='upper center',
                      bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True, shadow=True)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        fig.colorbar(self.scalarcm)

        if self.grid_off == True:
            ax.grid('off')
        if self.axis_off == True:
            ax.axis('off')

        return plt
        #plt.show()

    def plot_wulff_color(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        wulff_pt_list = self.wulff_pt_list
        on_wulff = self.on_wulff

        for plane in on_wulff:
            plane_color = plane[3]
            print plane_color, plane[0]
            pts = []
            for vertices in plane[2]:
                i = vertices[0]
                j = vertices[1]
                k = vertices[2]

                pts_vertices = [wulff_pt_list[i], wulff_pt_list[j], wulff_pt_list[k]]
                for n in range(1, 200):
                    pt1 = wulff_pt_list[i] + 0.005 * n * (wulff_pt_list[j]-wulff_pt_list[i])
                    pt2 = wulff_pt_list[i] + 0.005 * n * (wulff_pt_list[k]-wulff_pt_list[i])
                    pts_vertices.append(pt1)
                    pts_vertices.append(pt2)
                    for s in range(1, 200):
                        pt3 = (pt1 * s + pt2 * (200 - s)) * 0.005
                        pts_vertices.append(pt3)
                pts += pts_vertices

            ax.plot([x[0] for x in pts], [x[1] for x in pts], [x[2] for x in pts], color=plane_color, alpha=self.alpha)

            for line in plane[-1]:
                edge = [wulff_pt_list[line[0]], wulff_pt_list[line[1]]]
                ax.plot([x[0] for x in edge], [x[1] for x in edge], [x[2] for x in edge], 'k', lw=1)

        plt.gca().set_aspect('equal', adjustable='box')
        if self.show_area == True:
            ax.legend(self.color_proxy, self.miller_area, loc='upper left',
                      bbox_to_anchor=(-0.2, 1.05), fancybox=True, shadow=True)
        else:
            ax.legend(self.color_proxy, self.input_miller, loc='upper center',
                      bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True, shadow=True)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        fig.colorbar(self.scalarcm, alpha=self.alpha)

        # [normal, e_surf, normal_pt, dual_pt, color_plane, m_ind_orig, miller]
        e_range = max(self.bar_range) - min(self.bar_range)
        input_miller_ind = xrange(len(self.input_miller))
        if self.label_miller == True:
            for plane in self.normal_e_m:
                # normal pts on plane, add label there
                plot_m = 0
                unplot_miller = []
                for i in input_miller_ind:
                    if plane[-2] == i:
                        plot_m = 1
                    else:
                        unplot_miller.append(i)
                if plot_m == 0:
                    continue
                input_miller_ind = unplot_miller
                print 'zan~', plane[2], plane[4]
                normal_pts = [x * (plane[1] + 0.1 * e_range) for x in plane[0]]
                m_orig = '[' + str(plane[-1][0]) + str(plane[-1][1]) + str(plane[-1][2]) + ']'
                zdir = [plane[0][0], plane[0][1], plane[0][2]]
                ax.text(normal_pts[0], normal_pts[1], normal_pts[2], m_orig, zdir)
                if len(unplot_miller) == 0:
                    break

        if self.grid_off == True:
            ax.grid('off')
        if self.axis_off == True:
            ax.axis('off')

        plt.draw()

        return plt


    def plot_wulff_color_3v(self):
        fig1 = plt.figure(1)
        fig2 = plt.figure(2)
        fig3 = plt.figure(3)
        fig4 = plt.figure(4)
        #  top view
        ax1 = fig1.add_subplot(1, 1, 1, projection='3d', azim=0, elev=90)
        ax2 = fig2.add_subplot(1, 1, 1, projection='3d', azim=0, elev=90)
        ax3 = fig3.add_subplot(1, 1, 1, projection='3d', azim=0, elev=90)
        ax4 = fig4.add_subplot(1, 1, 1, projection='3d')
        wulff_pt_list = self.wulff_pt_list
        on_wulff = self.on_wulff

        for plane in on_wulff:
            plane_color = plane[3]
            print plane_color, plane[0]

            for vertices in plane[2]:
                i = vertices[0]
                j = vertices[1]
                k = vertices[2]
                pts_vertices = [wulff_pt_list[i], wulff_pt_list[j], wulff_pt_list[k]]
                v1 = pts_vertices[1] - pts_vertices[0]
                v2 = pts_vertices[2] - pts_vertices[0]

                data_test = array([wulff_pt_list[i], wulff_pt_list[j], wulff_pt_list[k]])
                Xs = data_test[:,0]
                Ys = data_test[:,1]
                Zs = data_test[:,2]
                # top view
                if abs(np.dot(np.cross(v1, v2),(0, 0, 1))) > 10e-10:
                    # print 'top'
                    ax1.plot_trisurf(Xs, Ys, Zs, color=plane_color, linewidth=0, alpha=self.alpha)

                # front view
                if abs(np.dot(np.cross(v1, v2),(1, 0, 0))) > 10e-10:
                    # print 'front'
                    ax2.plot_trisurf(Ys, Zs, Xs, color=plane_color, linewidth=0, alpha=self.alpha)

                # side view
                if abs(np.dot(np.cross(v1, v2),(0, 1, 0))) > 10e-10:
                    # print 'side'
                    ax3.plot_trisurf(Zs, Xs, Ys, color=plane_color, linewidth=0, alpha=self.alpha)

            for line in plane[-1]:
                edge = [wulff_pt_list[line[0]], wulff_pt_list[line[1]]]
                ax1.plot([x[0] for x in edge], [x[1] for x in edge], [x[2] for x in edge], 'k', lw=1.5)
                ax2.plot([x[1] for x in edge], [x[2] for x in edge], [x[0] for x in edge], 'k', lw=1.5)
                ax3.plot([x[2] for x in edge], [x[0] for x in edge], [x[1] for x in edge], 'k', lw=1.5)
                ax4.plot([x[0] for x in edge], [x[1] for x in edge], [x[2] for x in edge], 'k', lw=1.5)

        ax1.set_aspect('equal', adjustable='box')
        if self.show_area == True:
            ax1.legend(self.color_proxy, self.miller_area, loc='upper left',
                       bbox_to_anchor=(-0.2, 1.05), fancybox=True, shadow=True)
        else:
            ax1.legend(self.color_proxy, self.input_miller, loc='upper center',
                       bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True, shadow=True)
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_zlabel('z')
        ax1.set_title('top view', loc='right')
        fig1.colorbar(self.scalarcm, alpha=self.alpha)
        if self.grid_off == True:
            ax1.grid('off')
        if self.axis_off == True:
            ax1.axis('off')


        ax2.set_aspect('equal', adjustable='box')
        if self.show_area == True:
            ax2.legend(self.color_proxy, self.miller_area, loc='upper left',
                       bbox_to_anchor=(-0.2, 1.05), fancybox=True, shadow=True)
        else:
            ax2.legend(self.color_proxy, self.input_miller, loc='upper center',
                       bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True, shadow=True)
        ax2.set_xlabel('y')
        ax2.set_ylabel('z')
        ax2.set_zlabel('x')
        ax2.set_title('front view', loc='right')
        fig2.colorbar(self.scalarcm, alpha=self.alpha)
        if self.grid_off == True:
            ax2.grid('off')
        if self.axis_off == True:
            ax2.axis('off')

        ax3.set_aspect('equal', adjustable='box')
        if self.show_area == True:
            ax3.legend(self.color_proxy, self.miller_area, loc='upper left',
                       bbox_to_anchor=(-0.2, 1.05), fancybox=True, shadow=True)
        else:
            ax3.legend(self.color_proxy, self.input_miller, loc='upper center',
                       bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True, shadow=True)
        ax3.set_xlabel('z')
        ax3.set_ylabel('x')
        ax3.set_zlabel('y')
        ax3.set_title('side view', loc='right')
        fig3.colorbar(self.scalarcm, alpha=self.alpha)
        if self.grid_off == True:
            ax3.grid('off')
        if self.axis_off == True:
            ax3.axis('off')

        ax4.set_aspect('equal', adjustable='box')
        if self.show_area == True:
            ax4.legend(self.color_proxy, self.miller_area, loc='upper left',
                       bbox_to_anchor=(-0.2, 1.05), fancybox=True, shadow=True)
        else:
            ax4.legend(self.color_proxy, self.input_miller, loc='upper center',
                       bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True, shadow=True)
        ax4.set_xlabel('x')
        ax4.set_ylabel('y')
        ax4.set_zlabel('z')
        ax4.set_title('3d line view', loc='right')
        fig4.colorbar(self.scalarcm, alpha=self.alpha)
        if self.grid_off == True:
            ax4.grid('off')
        if self.axis_off == True:
            ax4.axis('off')

        plt.draw()
        return plt

    # with brewer color
    def plot_wulff_color_bw(self, azim=30, elev=60):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d', azim=azim, elev=elev)
        wulff_pt_list = self.wulff_pt_list
        on_wulff = self.on_wulff
        bw_color = self.brewer_color

        for plane in on_wulff:
            # get color from the order of original brewer number
            plane_color = bw_color[plane[4]]
            print plane_color, plane[0]
            pts = []
            for vertices in plane[2]:
                i = vertices[0]
                j = vertices[1]
                k = vertices[2]

                pts_vertices = [wulff_pt_list[i], wulff_pt_list[j], wulff_pt_list[k]]
                for n in range(1, 200):
                    pt1 = wulff_pt_list[i] + 0.005 * n * (wulff_pt_list[j] - wulff_pt_list[i])
                    pt2 = wulff_pt_list[i] + 0.005 * n * (wulff_pt_list[k] - wulff_pt_list[i])
                    pts_vertices.append(pt1)
                    pts_vertices.append(pt2)
                    for s in range(1, 200):
                        pt3 = (pt1 * s + pt2 * (200 - s)) * 0.005
                        pts_vertices.append(pt3)
                pts += pts_vertices

            ax.plot([x[0] for x in pts], [x[1] for x in pts], [x[2] for x in pts], color=plane_color, alpha=self.alpha)

            for line in plane[-1]:
                edge = [wulff_pt_list[line[0]], wulff_pt_list[line[1]]]
                ax.plot([x[0] for x in edge], [x[1] for x in edge], [x[2] for x in edge], 'k', lw=1)

        plt.gca().set_aspect('equal', adjustable='box')
        color_proxy = [plt.Rectangle((2, 2), 1, 1, fc=x, alpha=self.alpha) for x in self.brewer_color]
        if self.show_area == True:
            ax.legend(color_proxy, self.miller_area, loc='upper left',
                      bbox_to_anchor=(-0.2, 1.05), fancybox=True, shadow=True)
        else:
            ax.legend(color_proxy, self.input_miller, loc='upper center',
                      bbox_to_anchor=(0.5, 1.05), ncol=2, fancybox=True, shadow=True)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        ax1 = fig.add_axes([0.75, 0.15, 0.05, 0.70])
        cmap = colors.ListedColormap(self.brewer_color)
        cmap.set_over('0.75')
        cmap.set_under('0.25')
        bounds = self.input_esurf
        bounds.append(2 * bounds[-1] - bounds[-2])
        norm = colors.BoundaryNorm(bounds, cmap.N)
        cbar = colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm,
                                     boundaries=[0] + bounds + [10],
                                     extend='both',
                                     ticks=bounds,  # optional
                                     spacing='proportional',
                                     orientation='vertical')
        cbar.set_label('Surface Energies ($J/m^2$)')

        # [normal, e_surf, normal_pt, dual_pt, color_plane, m_ind_orig, miller]
        e_range = max(self.bar_range) - min(self.bar_range)
        input_miller_ind = xrange(len(self.input_miller))
        if self.label_miller == True:
            for plane in self.normal_e_m:
                # normal pts on plane, add label there
                plot_m = 0
                unplot_miller = []
                for i in input_miller_ind:
                    if plane[-2] == i:
                        plot_m = 1
                    else:
                        unplot_miller.append(i)
                if plot_m == 0:
                    continue
                input_miller_ind = unplot_miller
                print 'zan~', plane[2], plane[4]
                normal_pts = [x * (plane[1] + 0.1 * e_range) for x in plane[0]]
                m_orig = '[' + str(plane[-1][0]) + str(plane[-1][1]) + str(plane[-1][2]) + ']'
                zdir = [plane[0][0], plane[0][1], plane[0][2]]
                ax.text(normal_pts[0], normal_pts[1], normal_pts[2], m_orig, zdir)
                if len(unplot_miller) == 0:
                    break

        if self.grid_off == True:
            ax.grid('off')
        if self.axis_off == True:
            ax.axis('off')

        plt.draw()

        return plt
