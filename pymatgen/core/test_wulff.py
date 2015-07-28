__author__ = 'vivid0036'

from pymatgen import Structure, Lattice

from pymatgen.core.wulff_dual import wulff_3d
from pymatgen.io.smartio import CifParser

mo = Structure(Lattice.cubic(3.16), ["Mo", "Mo"],
                            [[0, 0, 0], [0.5, 0.5, 0.5]])
# li7 = Structure.from_file(get_path("relaxed_li7p3s11.cif"), primitive=False)
# e_surf_li7 = [0.21, 0.23, 0.25, 0.26, 0.34]
# miller_li7 = [[1, 0, 0], [1, 0, -1], [1, 0, 1], [0, 1, 0], [1, 1, -2]]
# wulff_li7 = wulff_3d(li7, miller_li7, e_surf_li7)
# wulff_li7.plot_wulff_color()
# plt_li7 = wulff_li7.plot_wulff_line()
# plt_li7.show()

e_surf_list = [3.34, 2.92, 3.24]

miller_list = [[0, 0, 1], [1, 1, 0], [1, 1, 1]]

wulff_mo = wulff_3d(mo, miller_list, e_surf_list)
# wulff_cu.plot_wulff_pts()

plt1 = wulff_mo.plot_wulff_line()


wulff_mo.plot_wulff_color()
plt1.show()