# -*- coding: utf-8 -*-
"""
parametric_image.py
One Ring to Bind them
- support notebook interface for BUFA graphic image

Started Thursday February 8, 2018
@author:
lanier4@illinois.edu
mradmstr514226508@gmail.com

Restarted - Redefined Sunday Seis de Mayo, 2018
"""
import warnings
warnings.filterwarnings('ignore')

import os
import sys
import numpy as np
import time

from IPython.display import display

sys.path.insert(1, '../src');
import z_plane as zp
import graphic_utility as gu;
import itergataters as ig
import numcolorpy as ncp


class BUFA:
    """ Behaviorally Unecpected Functional Algorithms BUFA

        Save, Load, view-set parameters display with Save AS output and Write Image AS features
            - notebook interface with dropdowns

        Wrappings for an Image:
            - Function definition with parameters
            - pixel numerical products definition - Three tightly coupled domains:
                - CGS, MKS or :) English defintion of Size and Resolution
                - class ComplexPlane definition
                    _center_point = CP
                    _zoom_factor = max(ZM, 1e-15)
                    _theta = theta
                    _n_rows = max(round(h), 1)
                    _n_cols = max(round(w), 1)
                - class ETA_iterator definition
                    _it_max = 32
                    _max_d = 12
            - pixel Color definition (partially implemented)
                - HSV = f(ET, theta, d)
                - RGB = f(ET | theta | d)
                - hybrid f(HSV, + RGB)

        With Functionality:
            - Resize
            - Reiterate
            - Recolor
    """

    def __init__(self, plane=zp.ComplexPlane()):

        """ ParametricImage constructor """
        self._plane = plane
        self.rnd_lambda = zp.rnd_lambda()
        print('ParametricImage happy happy initialization under construction', self.rnd_lambda)

