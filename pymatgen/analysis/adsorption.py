# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, print_function, unicode_literals
from __future__ import absolute_import

"""
This module provides classes used to enumerate surface sites
and to find adsorption sites on slabs
"""

import numpy as np
from six.moves import range

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


