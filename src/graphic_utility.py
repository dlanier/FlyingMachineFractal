# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 14:12:17 2016

@author: del
lanier4@illinois.edu
mradmstr514226508@gmail.com

"""
import time

import numpy as np
from PIL import Image as IP
from PIL import ImageColor as IC
import colorsys

def mat2graphic(Z):
    """ M, nClrs = mat2graphic(Z)
        Use all the transformation tricks to prepare input matrix Z
        for conversion to a viewable image.
    Args:
        Z:      real or complex (rows x xcols x 1) matrix
    Returns:
        M:      real (rows x xcols x 1) matrix (0 <= M <= 1)
    """
    M, nClrs = flat_index(np.abs(Z))
    return graphic_norm(M)
    
    
def graphic_norm(Z):
    """ rescale matrix z to distance (float) s.t.   
        0 <= z <= 1  (will include 0,1 if it has more than 1 value)
  
    Args:
        Z: is a real or complex two dimensional matrix
    
    Returns:
        Z: same size, real valued matrix with smallest member = 0, largest = 1
    """
    EPSILON = 1e-15
    I = np.abs(Z)
    I = I - I.min()
    
    return I / max(EPSILON, I.max())


def flat_index(float_mat):
    """ convert the input matrix to integers from 0 to number of unique values.
    
    Args:
        float_mat: two dimensional matrix.
        
    Return:
        float_mat: re-enumerated so that the matrix values are all sequential ints.
        n_colors:  number of unique values in the input / output matrix
    """
    rows = float_mat.shape[0]
    cols = float_mat.shape[1]

    float_mat = np.reshape(float_mat, (1, float_mat.size))
    ixA = np.argsort(float_mat)[0]
    
    current_value = float_mat[0, ixA[0]]
    enumeration_value = 0
    for ix in ixA:
        if float_mat[0,ix] != current_value:
            current_value = float_mat[0,ix]
            enumeration_value += 1
        float_mat[0,ix] = enumeration_value
    
    float_mat = np.array(np.reshape(float_mat, (rows, cols)))
    float_mat = np.int_(float_mat)
    n_colors = enumeration_value + 1
    
    return float_mat, n_colors


def resize_color_map(mp0, n_colors):
    """ givin a RGB colormap input return the same color order with n_colors number of colors
    """
    mp = np.zeros((n_colors,3))
    n_colors0 = mp0.shape[0]
    if n_colors0 != n_colors:
        tc = n_colors0 * n_colors
        x = np.linspace(1,tc, n_colors0)
        xq = np.linspace(1,tc, n_colors)
        mp[:,0] = np.interp(xq, x, mp0[:,0])
        mp[:,1] = np.interp(xq, x, mp0[:,1])
        mp[:,2] = np.interp(xq, x, mp0[:,2])
    return mp

def get_colormap_heat(n_colors):
    """ get at heat colormap
    """
    mp0 = resize_color_map(
        np.array([[1.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.0, 1.0]]),
        n_colors)
    
    return mp0

def int_mat_to_heat_image(int_arr, n_colors):
    """ convert a matrix of integers to an image with heat colormap
        where n_colors is the largest integer.
    """
    nRows = int_arr.shape[0]
    nCols = int_arr.shape[1]
    
    M = np.zeros((nRows, nCols, 3))
    mp = get_colormap_heat(n_colors+1)
    
    for r in range(0,nRows):
        for c in range(0,nCols):
            M[r, c, 0] = mp[int_arr[r, c], 0]
            M[r, c, 1] = mp[int_arr[r, c], 1]
            M[r, c, 2] = mp[int_arr[r, c], 2]
            
    return M


def hsv_2_rgb(h, s, v):
    """ http://stackoverflow.com/questions/24852345/hsv-to-rgb-color-conversion """
    if s == 0.0: return [v, v, v]
    i = int(h*6.) # XXX assume int() truncates!
    f = (h*6.)-i; p,q,t = v*(1.-s), v*(1.-s*f), v*(1.-s*(1.-f)); i%=6
    if i == 0: return [v, t, p]
    if i == 1: return [q, v, p]
    if i == 2: return [p, v, t]
    if i == 3: return [p, q, v]
    if i == 4: return [t, p, v]
    if i == 5: return [v, p, q]


def hsv_to_rgb(hsv_mat, rows, cols):
    """ rgb_mat = hsv_to_rgb(hsv_mat)
        using colorsys primatives.
    Args:
        hsv_mat: rows x cols x 3 hue-saturation-values image
    Returns:
        rgb_mat: rows x cols x 3 red-green-blue image
    """
    #rows = hsv_mat.shape[1]
    #cols = hsv_mat.shape[2]
    rgb_mat = np.zeros((3, rows, cols))
    for r in range(0, rows):
        for c in range(0, cols):
            h = hsv_mat[0, r, c]
            s = hsv_mat[1, r, c]
            v = hsv_mat[2, r, c]
            red, green, blue = colorsys.hsv_to_rgb(h, s, v)
            rgb_mat[0, r, c] = red
            rgb_mat[1, r, c] = green
            rgb_mat[2, r, c] = blue
    #rgb_mat = IP.fromarray(rgb_mat, mode='RGB')
    return rgb_mat

def mats3_as_hsv_2rgb(H, S, V):
    """ I = mats3_as_hsv_2rgb(H, S, V)
        combine three normaized matrices into an hsv image
        'normalized' means     0 <= M <= 1   where M is H, S, or V
    Args:
        H:  rows x cols x 1 normalized matrix
        S:  rows x cols x 1 normalized matrix
        V:  rows x cols x 1 normalized matrix
    Returns:
        I:  rows x cols x 3  hue-saturation-values image
    """
    n_rows = H.shape[0]
    n_cols = H.shape[1]
    I = IP.new('RGB', (n_cols, n_rows))

    for row in range(0, I.height):
        for col in range(0, I.width):
            #P = tuple(np.int_([H[row, col] * 255, S[row, col] * 255, V[row, col] * 255]))
            red, green, blue = colorsys.hsv_to_rgb(H[row, col], S[row, col], V[row, col])
            red = int(np.round(red*255))
            green = int(np.round(green*255))
            blue = int(np.round(blue*255))
            P = (red, green, blue)
            I.putpixel((col, row), P)
            
    return I

def paint_ET_seq(ET, max_w=800):
    n_max = ET.max()
    ET_seq = np.int_(ET == 0)
    for n in range(1, n_max):
        ET_seq = np.concatenate((ET_seq, np.int_(ET==n) ))

    return ET_seq

def Z_ET_to_show(Z, ET):
    """ quick showable graphic of escape-time calculation results """
    V = 1 - graphic_norm(Z)
    S = 1 - graphic_norm(ET)
    H = graphic_norm(np.arctan2(np.imag(Z), np.real(Z)))
    return mats3_as_hsv_2rgb(H, S, V)


def Z_ET_to_show_II(Z, ET):
    """ quick showable graphic of escape-time calculation results """
    H = graphic_norm(Z)
    S = graphic_norm(ET)
    V = 1 - S
    I = mats3_as_hsv_2rgb(H, S, V)
    return I


def ET_as_grayscale(ET):
    """ I = ET_as_grayscale(ET)
    Args:
        ET: rows x cols x 1 numerical (real) matrix
    Returns:
        I:  rows x cols x 3 image
    """
    Iet = IP.fromarray(graphic_norm(ET) * 255)
    return Iet

def mat_to_gray(V):
    """ I = mat_to_gray(V)
        matrix of values V, converted to a gray scale image
    Args:
        V:  rows x cols x 1 numerical matrix
    Returns:
        I:  rows x cols x 3 grayscale image
    """
    V = mat2graphic(V)
    n_rows = V.shape[0]
    n_cols = V.shape[1]
    I = IP.new('RGB', (n_cols, n_rows))
    for row in range(0, I.height):
        for col in range(0, I.width):
            P = tuple(np.int_(
                [V[row, col] * 255, V[row, col] * 255, V[row, col] * 255]))
            I.putpixel((col, row), P)
    return I


def mat_to_red(V):
    R_max = 255
    R_floor = 180
    G_max = 250
    G_floor = 30
    B_max = 250
    B_floor = 30
    return mat_to_RGB(V, R_max, G_max, B_max, R_floor, G_floor, B_floor)


def mat_to_green(V):
    R_max = 250
    R_floor = 30
    G_max = 255
    G_floor = 130
    B_max = 250
    B_floor = 30
    return mat_to_RGB(V, R_max, G_max, B_max, R_floor, G_floor, B_floor)


def mat_to_blue(V):
    R_max = 250
    R_floor = 30
    G_max = 250
    G_floor = 30
    B_max = 255
    B_floor = 130
    return mat_to_RGB(V, R_max, G_max, B_max, R_floor, G_floor, B_floor)


def mat_to_RGB(V, R_max, G_max, B_max, R_floor=0, G_floor=0, B_floor=0):
    """ I = mat_to_gray(V)
        matrix of values V, converted to a gray scale image
    Args:
        V:  rows x cols x 1 numerical matrix
    Returns:
        I:  rows x cols x 3 grayscale image
    """
    R = R_max - R_floor
    G = G_max - G_floor
    B = B_max - B_floor
    V = 1 - mat2graphic(V)
    n_rows = V.shape[0]
    n_cols = V.shape[1]
    I = IP.new('RGB', (n_cols, n_rows))
    for row in range(0, I.height):
        for col in range(0, I.width):
            P = tuple(np.int_(
                [R_floor + V[row, col] * R, G_floor + V[row, col] * G, B_floor + V[row, col] * B]))
                
            I.putpixel((col, row), P)
    return I

def get_a_snowscreen(h, w):
    """ color_snow_image = get_a_snowscreen(h, w)
        get a pseudo-random rgb PIL image
    Args:       h, w == height and width of image
    Returns:    color_snow h x w x 3 image
    """
    snow_screen = IP.new('RGB', (w, h))
    for row in range(0, snow_screen.height):
        for col in range(0, snow_screen.width):
            snow_screen.putpixel((col, row), tuple(np.int_(np.random.random(3)*255)))
        
    return snow_screen

def show_matrix(A, n_dec=3):
    """ Display a matrix (nicely) with a fixed number of decimal places.
    Args:
    """
    format_str = '%' + '0.0%df\t'%(n_dec)
    n_rows = A.shape[0]
    n_cols = 1
    if A.size > sum(A.shape): n_cols = A.shape[1]
    for row in range(0, n_rows):
        s = ''
        for col in range(0, n_cols):
            s = s + format_str%(A[row,col])
        print(s)
    return

def plane_gradient(X):
    """ DX, DY = plane_gradient(X) 
    Args:
        X:          matrix
    Returns:
        DX:         gradient in X direction
        DY:         gradient in Y direction
    """
    n_rows = X.shape[0]
    n_cols = X.shape[1]
    DX = np.zeros(X.shape)
    DY = np.zeros(X.shape)
    for r in range(0, n_rows):
        xr = X[r, :]
        for c in range(0, n_cols - 1):
            DX[r,c] = xr[c+1] - xr[c]
            
    for c in range(0, n_cols):
        xc = X[:, c]
        for r in range(0, n_rows -1):
            DY[r, c] = xc[r+1] - xc[r]
        
    return DX, DY

def grad_Im(X):
    """
    Args:
        X:               matrix
    Returns:
        Gradient_Image:  positive matrix representation of the X-Y gradient of X
    """
    DX, DY = plane_gradient(X)
    return graphic_norm(DX + DY * 1j)


def integer_to_alphabet_number(N):
    """ convert integer to spreadsheet column lettering 
    
    Args:    an integer representation of the letters like 136964
    
    Returns: a string of capital letters like 'GTOV'
    """
    B = 26
    ALPH = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    max_mag = 100
    M = 0

    while (N > B**M) & (M <= max_mag):
        M += 1
        
    S = []
    for nada in range(0, M):
        M -= 1
        D = np.int_(np.floor(N / B**M))
        N = N - D * B**M
        S.append(ALPH[D-1])
    
    alpha_N = ''.join(S)
    
    return alpha_N

def alphabet_number_to_integer(N):
    """ convert spreadsheet column lettering to an integer 
    
    Args:    a string of capital letters like 'GTOV'
    
    Returns: an integer representation of the letters like 136964
    """
    ALPH = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    alphanumeric_dictionary = dict(zip(ALPH, np.arange(1, len(ALPH)+1)))
    N = list(N)[::-1]
    numeric_N = 0
    for p in range(0, len(N)):
        numeric_N += alphanumeric_dictionary[N[p]] * 26**p
    return numeric_N

def uniq_time_id_string(time_units=1e6, Warn_The_World=True):
    """ time in (1 <= time_units <= 1e9) as integer to string...  """
    if Warn_The_World & (time_units > 1e9): print('uniq_time_id_string:\t You Fool!')
    N = int(time.time() * np.maximum(np.minimum(time_units, 1e9), 1))
    uniq_id = integer_to_alphabet_number(N)
    return uniq_id, N

def now_name(n_digits=9, prefi_str=None, suffi_str=None):
    """ get a human readable time stamp name unique up to 1 / 1e9 seconds """
    t0 = time.time()
    t_dec = t0 - np.floor(t0)
    ahora_nombre = time.strftime("%a_%d_%b_%Y_%H_%M_%S",
                                 time.localtime()) + '_' + ('%0.09f'%(t_dec))[2:n_digits+2]
    if prefi_str is None: prefi_str = ''
    if suffi_str is None: suffi_str = ''
    return prefi_str + ahora_nombre + suffi_str

class hsv_object:

    def __init__(self, Z0, Z, ET):
        self._Z0 = Z0
        self._Z = Z
        self._ET = ET
        self._D = np.abs(Z-Z0)
        self._R = np.arctan2(np.imag(Z-Z0), np.real(Z-Z0))

        self.V = 1 - mat2graphic(self._D)
        self.S = 1 - mat2graphic(self._ET)
        self.H = mat2graphic(self._R)

    def get_image(self):

        # self.S = 1 - mat2graphic(self._D)
        # self.V = 1 - mat2graphic(self._ET)
        # self.H = mat2graphic(self._R)

        return mats3_as_hsv_2rgb(self.H, self.S, self.V)

# def Z_ET_to_show(Z, ET):
#     """ quick showable graphic of escape-time calculation results """
#     V = 1 - graphic_norm(Z)
#     S = 1 - graphic_norm(ET)
#     H = graphic_norm(np.arctan2(np.imag(Z), np.real(Z)))
#     return mats3_as_hsv_2rgb(H, S, V)


