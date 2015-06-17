__author__ = 'richard'

import os
from pymongo import MongoClient
from fireworks.core.firework import FireWork, Workflow
from fireworks.core.launchpad import LaunchPad
from pymatgen.io.vaspio_set import MPVaspInputSet
from pymatgen.transformations.standard_transformations import RemoveSpeciesTransformation
from pymacy.fireworks_scratch.anode_task_zix009 import RunCustodianTask, VaspDBInsertTask
from pymacy.fireworks_scratch.anode_task_zix009 import VASPInputFromID, get_id_list_mo