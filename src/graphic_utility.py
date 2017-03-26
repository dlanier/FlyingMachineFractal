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
    """ M = mat2graphic(Z)
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

def hsv_to_rgb(hsv_mat):
    """ rgb_mat = hsv_to_rgb(hsv_mat)
        using colorsys primatives.
    Args:
        hsv_mat: rows x cols x 3 hue-saturation-values image
    Returns:
        rgb_mat: rows x cols x 3 red-green-blue image
    """
    rows = hsv_mat.shape[1]
    cols = hsv_mat.shape[2]
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
    # I = IP.new('RGB', (n_cols, n_rows))
    I = IP.new('HSV', (n_cols, n_rows))
    for row in range(0, I.height):
        for col in range(0, I.width):
            P = tuple(np.int_(
                [H[row, col] * 255, S[row, col] * 255, V[row, col] * 255]))
            I.putpixel((col, row), P)
    return hsv_to_rgb(I)

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

def Z_ET_to_show(Z, ET):
    """ quick showable graphic of escape-time calculation results """
    H = 1 - graphic_norm(Z)
    V = 1 - graphic_norm(ET)
    S = graphic_norm(np.arctan2(np.imag(Z), np.real(Z)))
    I = mats3_as_hsv_2rgb(H, S, V)
    return I

def Z_ET_to_show_II(Z, ET):
    """ quick showable graphic of escape-time calculation results """
    H = graphic_norm(Z)
    S = graphic_norm(ET)
    V = 1 - S
    I = mats3_as_hsv_2rgb(H, S, V)
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
    ahora_nombre = time.strftime("%a_%d_%b_%Y_%H_%M_%S", time.localtime()) + '_' + ('%0.09f'%(t_dec))[2:n_digits+2]
    if prefi_str is None: prefi_str = ''
    if suffi_str is None: suffi_str = ''
    return prefi_str + ahora_nombre + suffi_str
