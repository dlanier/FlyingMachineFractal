"""
complex plane functions and class for getting pixels in three different measuring systems
"""
import sys
import numpy as np
from zexhibit import get_aligned_dict_string

Zomm_factor_min = 1000 * sys.float_info.min

_default_dict = {'center_point': 0.0 + 0.0j,
                'zoom_factor': 0.5,
                'theta': 0.0,
                'n_rows':400,
                'n_cols':400}

def get_frame_from_dict(def_dict):
    """ complex_frame, def_dict = get_frame_from_dict(def_dict)
        legacy wrapper function.
    Args:
        def_dict: definition dictionary with keys:
                    'center_point', 'zoom', 'theta', 'n_rows', 'n_cols'
    Returns:
        complex_frame:
        def_dict:
    """
    complex_frame = get_complex_frame(
        def_dict['center_point'],
        def_dict['zoom_factor'],
        def_dict['theta'],
        def_dict['n_rows'],
        def_dict['n_cols'])

    return complex_frame, def_dict

def get_complex_frame(CP, ZM, theta, h=1, w=1):
    """ get the complex numbers at ends and centers of a frame  """
    frame_dict = {'center_point':CP}
    if w >= h:
        frame_dict['top_center'] = np.exp(1j*(np.pi/2 + theta))/ZM
        frame_dict['right_center'] = (w/h) * np.exp(1j * theta) / ZM
    else:
        frame_dict['top_center'] = (h/w) * np.exp(1j*(np.pi/2 + theta)) / ZM
        frame_dict['right_center'] = np.exp(1j * theta) / ZM

    frame_dict['bottom_center'] = frame_dict['top_center'] * -1
    frame_dict['left_center'] = frame_dict['right_center'] * -1
    frame_dict['upper_right'] = frame_dict['right_center'] + frame_dict['top_center']
    frame_dict['bottom_right'] = frame_dict['right_center'] + frame_dict['bottom_center']

    frame_dict['upper_left'] = frame_dict['left_center'] + frame_dict['top_center']
    frame_dict['bottom_left'] = frame_dict['left_center'] + frame_dict['bottom_center']

    for k in frame_dict.keys():
        frame_dict[k] = frame_dict[k] + CP

    return frame_dict

class ComplexPlain:
    """                         parameterized grid of complex numbers
    Args:
        CP:                     self._center_point    -- complex vector from origin to center of grid
        ZM:                     self._zoom_factor     -- Magnify (Zoom IN as ZM increases)
        theta:                  self._theta           -- Counter Clockwise rotation of the plane
        h:                      self._n_rows          -- number of rows in the grid
        w:                      self._n_cols          -- number of columns in the grid

    methods:
        display_self:           command line printout of the self definition parameters (Args)
        get_complex_axes:       grid center arrays of complex vectors
        get_complex_col:        column array of complex vectors
        get_complex_row:        row array of complex vectors
        get_complex_pixels:     matrix of complex numbers == the grid
        get_escape_bound:       get an escpe distance based on the grids corner to corner vector length
        get_parameters_dict:    get the self definition parameters (Args) as a python dict
        get_rails:              top and bottom arrays of complex vectors
        get_styles:             left and right arrays of complex vectors
        load_dict:              re-initialize the object with new set of definition parameters (Args)
        , CP=0.0+0.0*1j, ZM=1.0, theta=0.0, h=5, w=5
    """

    def __init__(self, plain_parameters=_default_dict):
        self._center_point = plain_parameters['center_point']
        self._zoom_factor = max(plain_parameters['zoom_factor'], Zomm_factor_min)
        self._theta = plain_parameters['']
        self._n_rows = max(round(plain_parameters['n_rows']), 1)
        self._n_cols = max(round(plain_parameters['n_cols']), 1)

    def display_self(self):
        pd = self.get_parameters_dict()
        s = get_aligned_dict_string(pd)
        print(s)
    
    def load_dict(self, parameters_dict):
        """ self.load_dict(parameters_dict) """
        if 'center_point' in parameters_dict:
            self._center_point = parameters_dict['center_point']
        if 'zoom_factor' in parameters_dict:
            self._zoom_factor = parameters_dict['zoom_factor']
        if 'theta' in parameters_dict:
            self._theta = parameters_dict['theta']
        if 'n_rows' in parameters_dict:
            self._n_rows = parameters_dict['n_rows']
        if 'n_cols' in parameters_dict:
            self._n_cols = parameters_dict['n_cols']

    def get_parameters_dict(self):
        """ parameters_dict = self.get_parameters_dict() """
        parameters_dict = {}
        parameters_dict['center_point'] = self._center_point
        parameters_dict['zoom_factor'] = self._zoom_factor
        parameters_dict['theta'] = self._theta
        parameters_dict['n_rows'] = self._n_rows
        parameters_dict['n_cols'] = self._n_cols    
        return parameters_dict

    def get_escape_bound(self, boundry_scale=12):
        """ escape time algorithm best infinity safe iteration distance """
        corner_scale = max(self._n_rows/self._n_cols, self._n_cols/self._n_rows)
        return corner_scale * boundry_scale / self._zoom_factor

    def get_complex_axes(self):
        """ horiz_axis, vert_axis = self.get_complex_axes() """
        frame_dict = get_complex_frame(self._center_point,
                                       self._zoom_factor, self._theta, self._n_rows, self._n_cols)
        vert_axis = np.linspace(frame_dict['top_center'],
                                frame_dict['bottom_center'], self._n_cols) + 0.0j
        horiz_axis = np.linspace(frame_dict['left_center'],
                                 frame_dict['right_center'], self._n_cols) + 0.0j

        return horiz_axis, vert_axis

    def get_rails(self):
        """ top_rail, bottom_rail = self.get_styles() """
        frame_dict = get_complex_frame(self._center_point,
                                       self._zoom_factor, self._theta, self._n_rows, self._n_cols)
        top_rail = np.linspace(frame_dict['upper_left'],
                               frame_dict['upper_right'], self._n_cols) + 0.0j
        bottom_rail = np.linspace(frame_dict['bottom_left'],
                                  frame_dict['bottom_right'], self._n_cols) + 0.0j

        return top_rail, bottom_rail

    def get_styles(self):
        """ left_style, right_style = self.get_styles() """
        frame_dict = get_complex_frame(self._center_point,
                                       self._zoom_factor, self._theta, self._n_rows, self._n_cols)
        left_style = np.linspace(frame_dict['upper_left'],
                                 frame_dict['bottom_left'], self._n_rows) + 0.0j
        right_style = np.linspace(frame_dict['upper_right'],
                                  frame_dict['bottom_right'], self._n_rows) + 0.0j

        return left_style, right_style

    def get_complex_row(self, row_number):
        """ row_vectors = self.get_complex_row() """
        left_style, right_style = self.get_styles()

        return np.linspace(left_style[row_number], right_style[row_number], self._n_cols) + 0.0j

    def get_complex_col(self, col_number):
        """ col_vectors = self.get_complex_col() """
        top_rail, bottom_rail = self.get_styles()

        return np.linspace(top_rail[col_number], bottom_rail[col_number], self._n_rows) + 0.0j

    def get_complex_pixels(self):
        """ complex_pixels = self.get_complex_pixels() """
        left_style, right_style = self.get_styles()
        complex_pixels = np.zeros((self._n_rows,
                                   self._n_cols)) + np.zeros((self._n_rows, self._n_cols)) * 1j

        for k in range(0, self._n_rows):
            complex_pixels[k, :] = np.linspace(left_style[k], right_style[k], self._n_cols)

        return complex_pixels

""" end def class ComplexPlain: """
