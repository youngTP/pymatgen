## for Surface Energy Calculation
from __future__ import division, unicode_literals

__author__ = "Zihan XU"
__version__ = "0.1"
__email__ = "vivid0036@gmail.com"
__date__ = "6/2/15"

import itertools
import json
import os
import os.path as pth

from pymongo import MongoClient
from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.utilities import fw_utilities
from fireworks import explicit_serialize
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.io.vaspio_set import MPVaspInputSet, DictVaspInputSet
from custodian.custodian import Custodian
from custodian.vasp.jobs import VaspJob
from matgendb.creator import VaspToDbTaskDrone

from pymatgen.core.surface import Slab, SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.matproj.rest import MPRester


"""
for VASP Input
"""
class MPSlabVaspInputSet(DictVaspInputSet):
    """
    Class for writing a slab vasp run with LDA method.

    Args:
        user_incar_settings(dict): A dict specifying additional incar
            settings, default to None
        kpoints0: specify kpts[0] = kpoints0, default to []
        k_product: kpts[0][0]*a. Decide k density without kpoint0,
        default to 45
        potcar_functional: default to LDA
        bulk: default to False

        **kwargs:
            Other kwargs supported by :class:`DictVaspInputSet`.
    """
    def __init__(self, user_incar_settings=None, kpoints0=[],
                  k_product=45, potcar_functional='LDA', bulk=False, **kwargs):
        vis = MPVaspInputSet().as_dict()
        DictVaspInputSet.__init__(self, "MaterialsProject Slab", vis["config_dict"],
                                  **kwargs)
        if bulk:
             self.incar_settings.update(
                 {"NPAR": 4})
        else:
            self.incar_settings.update(
                {"AMIN": 0.01, "AMIX": 0.2, "BMIX": 0.001, "ISIF": 2,
                 "NPAR": 4, "NELMIN": 8})
        self.user_incar_settings = user_incar_settings or {}
        if user_incar_settings:
            self.incar_settings.update(user_incar_settings)

        self.k_product = k_product
        self.kpoints0 = kpoints0
        self.potcar_functional = potcar_functional
        self.bulk = bulk


    def get_kpoints(self, structure):
        """
        kpoint0 is the first consideration,
        k_product is a second choice, default to 40
        """
        kpt = super(MPSlabVaspInputSet, self).get_kpoints(structure)
        kpt.comment = "Automatic mesh"
        kpt.style = 'Gamma'

        # use k_product to calculate kpoints, k_product = kpts[0][0] * a
        abc = structure.lattice.abc
        kpt_calc = [int(self.k_product/abc[0]+0.5),
                    int(self.k_product/abc[1]+0.5), 1]
        self.kpt_calc = kpt_calc
        # calculate kpts (c direction) for bulk. (for slab, set to 1)
        if self.bulk:
            kpt_calc[2] = int(self.k_product/abc[2]+0.5)

        # kpoint0 is prior to k_product
        if self.kpoints0:
            kpt.kpts[0] = self.kpoints0
            print "kpoint0", kpt.kpts
        else:
            kpt.kpts[0] = kpt_calc
            print 'kpt_calc:\n', kpt_calc

        return kpt

    def get_incar(self, structure):
        abc = structure.lattice.abc
        kpt_calc = [int(self.k_product/abc[0]+0.5),
                    int(self.k_product/abc[1]+0.5), int(self.k_product/abc[1]+0.5)]

        if self.kpoints0:
            kpts = self.kpoints0
        else:
            kpts = kpt_calc

        if kpts[0]<5 and kpts[1]<5:
            if not self.bulk:
                self.incar_settings.update(
                    {"ISMEAR": 0})
            else:
                if kpts[2]<5:
                    self.incar_settings.update(
                        {"ISMEAR": 0})
        if self.user_incar_settings:
                self.incar_settings.update(self.user_incar_settings)

        incr = super(MPSlabVaspInputSet, self).get_incar(structure)

        return incr

    def as_dict(self):
        d = super(MPSlabVaspInputSet, self).as_dict()
        d.update({
            "kpoints0": self.kpoints0,
            "potcar_functional": self.potcar_functional,
            "user_incar_settings": self.user_incar_settings
        })
        return d

    def from_dict(cls, d):
        return cls(kpoints0=d.get("kpoints0", []),
                   user_incar_settings=d.get("user_incar_settings", None),
                   potcar_functional=d.get("potcar_functional", None))

    def get_all_vasp_input(self, structure):
        """
        Returns all input files as a dict of {filename: vaspio object}

        Args:
            structure (Structure/IStructure): Structure to generate vasp
                input for.
        Returns:
            dict of {filename: file_as_string}, e.g., {'INCAR':'EDIFF=1e-4...'}
        """
        data = {'INCAR': self.get_incar(structure),
             'KPOINTS': self.get_kpoints(structure),
             'POSCAR': self.get_poscar(structure),
             'POTCAR': self.get_potcar(structure)}
        return data


"""
Firework tasks
"""

def get_input_mp(element, miller_index, api_key, min_slab_size=10, min_vacuum_size=10,
                 symprec=0.001, angle_tolerance=5):
    """
    element: str, element name of Metal
    miller_index: hkl, e.g. [1, 1, 0]
    api_key: to get access to MP DB
    """
    # This initializes the REST adaptor. Put your own API key in.
    # e.g. MPRester("QMt7nBdIioOVySW2")
    mprest = MPRester(api_key)
    #first is the lowest energy one
    prim_unit_cell = mprest.get_structures(element)[0]
    spa = SpacegroupAnalyzer(prim_unit_cell,  symprec=symprec,
                             angle_tolerance=angle_tolerance)
    conv_unit_cell = spa.get_conventional_standard_structure() 
    slab_metal = SlabGenerator(conv_unit_cell, miller_index, min_slab_size, min_vacuum_size, 
                               lll_reduce=False, center_slab=False, primitive=False)
    oriented_u_cell = slab_metal.oriented_unit_cell
    slab_cell = slab_metal.get_slab()



    return oriented_u_cell, slab_cell

@explicit_serialize
class WriteSurfVaspInput(FireTaskBase):
    """writes VASP inputs given elements, hkl,  """

    required_params = ["element", "miller_index", "api_key"]
    optional_params = ["min_slab_size", "min_vacuum_size",
                       "symprec", "angle_tolerance", "user_incar_settings",
                       "k_product","potcar_functional"]

    def run_task(self, fw_spec):
        dec = MontyDecoder()
        element = dec.process_decoded(self.get("element"))
        miller_index = dec.process_decoded(self.get("miller_index"))
        api_key = dec.process_decoded(self.get("api_key"))
        min_slab_size= dec.process_decoded(self.get("min_slab_size", 10))
        min_vacuum_size = dec.process_decoded(self.get("min_vacuum_size", 10))
        symprec = dec.process_decoded(self.get("symprec", 0.001))
        angle_tolerance = dec.process_decoded(self.get("angle_tolerance", 5))
        user_incar_settings = dec.process_decoded(self.get("user_incar_settings", 
                                                           {'ISIF': 2, 'EDIFFG':  -0.05,'EDIFF': 0.0001,
                                                            'ISMEAR': 1,'AMIX': 0.1,'BMIX': 0.0001, 
                                                            'AMIX_MAG': 0.4, 'BMIX_MAG': 0.0001, 
                                                            'NPAR':4, 'SIGMA': 0.05}))
        k_product = dec.process_decoded(self.get("k_product", 50))
        potcar_functional = dec.process_decoded(self.get("potcar_fuctional", 'LDA'))
        
        
        input_structures = get_input_mp(element, miller_index, api_key, min_slab_size, 
                                        min_vacuum_size,symprec, angle_tolerance)

        orient_u_cell = input_structures[0]
        slab_cell = input_structures[1]
        mplb_u = MPSlabVaspInputSet(user_incar_settings=user_incar_settings, k_product=k_product, 
                                    potcar_funtional=potcar_functional, bulk = True)
        mplb_u.write_input(orient_u_cell, '%s_ucell_k%s' %element %k_product)
        
        mplb_s = MPSlabVaspInputSet(user_incar_settings=user_incar_settings, k_product=k_product, 
                                    potcar_funtional=potcar_functional, bulk = False)
        mplb_s.write_input(slab_cell, '%s_scell_k%s' %element %k_product)
        
@explicit_serialize
class RunCustodianTask(FireTaskBase):
    """Runs Custodian."""

    required_params = ["jobs"]
    optional_params = ["custodian_params"]

    def run_task(self, fw_spec):

        fw_env = fw_spec.get("_fw_env", {})
        cust_params = self.get("custodian_params", {})
        if fw_env.get('scratch_root'):
            cust_params['scratch_dir'] = os.path.expandvars(
                fw_env['scratch_root'])

        dec = MontyDecoder()
        #handlers = dec.process_decoded(self['handlers'])
        jobs = dec.process_decoded(self['jobs'])
        #validators = [VasprunXMLValidator()]
        handlers = [VaspErrorHandler(), MeshSymmetryErrorHandler(),
                    UnconvergedErrorHandler(), NonConvergingErrorHandler(),
                    PotimErrorHandler()]

        c = Custodian(handlers, jobs, max_errors=10, **cust_params)
        output = c.run()

        return FWAction(stored_data=output)


@explicit_serialize
class VaspDBInsertTask(FireTaskBase):

    required_params = ["host", "port", "user", "password", "database",
                       "collection", "struct_type", "miller_index"]

    def run_task(self, fw_spec):

        dec = MontyDecoder()
        miller_index = dec.process_decoded(self.get("miller_index"))
        struct_type = dec.process_decoded(self.get("struct_type"))

        if not self["host"]:
            self["host"] = "127.0.0.1"

        if not self["port"]:
            self["port"] = 27017

        if not self["database"]:
            self["database"] = "vasp"

        if not self["collection"]:
            self["collection"] = "tasks"

        try:
            with open("FINAL_LOC") as f:
                loc = f.read()
        except IOError:
            loc = "."

        with open(pth.join("/projects/ong-group/repos/pymacy/resources/db/",
                            "db.json")) as f:
            d = json.load(f)
            drone = VaspToDbTaskDrone(host=self["host"], port=self["port"],
                                      user=self["user"], password=self["password"],
                                      database=self["database"], collection=self["collection"],
                                      mapi_key=d["mapi_key"],
                                      additional_fields={"author": os.environ.get("USER"),
                                                         "type": struct_type,
                                                         "miller index": miller_index},
                                      use_full_uri=False)
            drone.assimilate(loc)

