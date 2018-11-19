# -*- coding: utf-8 -*-
"""
numcolorpy.py
Created Saturday April 22 2017

@author: del
lanier4@illinois.edu
mradmstr514226508@gmail.com

import numcolorpy as ncp
"""
import time

import numpy as np
from PIL import Image as IP

from PIL import ImageColor as IC
import colorsys

def range_norm(Z, lo=0.0, hi=1.0):
    """ normaize input matrix Z within a lo - hi range 
    """
    I = graphic_norm(Z)
    hi = max(min(hi, 1.0), 0.0)
    lo = min(max(lo, 0.0), 1.0)
    low_fence = min(hi, lo)
    hi_fence = max(hi, lo)

    if low_fence == hi_fence:
        return I
    
    v_span = hi_fence - low_fence
    I = I * v_span + low_fence
    
    return I


def etg_norm(Z0, Z, ET):
    """ Zd, Zr, ETn = etg_norm(Z0, Z, ET); Graphically usable matrices from escape time algorithm result 
    """
    ETn = mat2graphic(ET)
    Zv = Z - Z0
    Zd = mat2graphic(Zv)
    Zr = mat2graphic(np.arctan2(np.imag(Zv), np.real(Zv)))
    return Zd, Zr, ETn


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


def gray_mat(V):
    n_rows = V.shape[0]
    n_cols = V.shape[1]
    V = V * 255
    I = IP.new('RGB', (n_cols, n_rows))
    for row in range(0, I.height):
        for col in range(0, I.width):
            P = tuple(np.int_([V[row, col], V[row, col], V[row, col]]))
            I.putpixel((col, row), P)
    return I


def rgb_2_hsv_mat(H, S, V):
    n_rows = H.shape[0]
    n_cols = H.shape[1]
    I = IP.new('RGB', (n_cols, n_rows))
    
    for row in range(0, I.height):
        for col in range(0, I.width):
            red, green, blue = colorsys.hsv_to_rgb(H[row, col], S[row, col], V[row, col])
            red = int(np.round( red * 255 ))
            green = int(np.round( green * 255 ))
            blue = int(np.round( blue * 255 ))
            P = (red, green, blue)
            I.putpixel((col, row), P)
            
    return I


def mat_to_gray(V, max_v=255, min_v=0):
    
    R_max = max(min(max_v, 255), 0)
    R_floor = min(max(min_v, 0), R_max)
    G_max = max(min(max_v, 255), 0)
    G_floor = min(max(min_v, 0), G_max)
    B_max = max(min(max_v, 255), 0)
    B_floor = min(max(min_v, 0), B_max)

    return mat_to_Shade(V, R_max, G_max, B_max, R_floor, G_floor, B_floor)


def mat_to_red(V):
    R_max = 255
    R_floor = 180
    G_max = 250
    G_floor = 30
    B_max = 250
    B_floor = 30
    return mat_to_Shade(V, R_max, G_max, B_max, R_floor, G_floor, B_floor)


def mat_to_green(V):
    R_max = 250
    R_floor = 30
    G_max = 255
    G_floor = 130
    B_max = 250
    B_floor = 30
    return mat_to_Shade(V, R_max, G_max, B_max, R_floor, G_floor, B_floor)


def mat_to_blue(V):
    R_max = 250
    R_floor = 30
    G_max = 250
    G_floor = 30
    B_max = 255
    B_floor = 130
    return mat_to_Shade(V, R_max, G_max, B_max, R_floor, G_floor, B_floor)


def mat_to_Shade(V, R_max, G_max, B_max, R_floor=0, G_floor=0, B_floor=0):
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
    V = graphic_norm(V)
    n_rows = V.shape[0]
    n_cols = V.shape[1]
    I = IP.new('RGB', (n_cols, n_rows))
    for row in range(0, I.height):
        for col in range(0, I.width):
            P = tuple(np.int_(
                [R_floor + V[row, col] * R, G_floor + V[row, col] * G, B_floor + V[row, col] * B]))
                
            I.putpixel((col, row), P)
    return I


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


def normat_hsv_intrgb(H, S, V, H_max=1.0, H_min=0.0, S_max=1.0, S_min=0.0, V_max=1.0, V_min=0.0):
    """ I = normat_hsv_intrgb(H, S, V, H_max=1.0, H_min=0.0, S_max=1.0, S_min=0.0, V_max=1.0, V_min=0.0)
        Three normaized matrices as hsv image converted to rgb
        'normalized' means     0 <= M <= 1   where M is H, S, or V
    Args:
        H:  rows x cols x 1 normalized matrix
        S:  rows x cols x 1 normalized matrix
        V:  rows x cols x 1 normalized matrix
    Returns:
        I:  rows x cols x 3  hue-saturation-values image
    """
    H_mul = H_max - H_min
    S_mul = S_max - S_min
    V_mul = V_max - V_min
    n_rows = H.shape[0]
    n_cols = H.shape[1]
    I = IP.new('RGB', (n_cols, n_rows))
    
    for row in range(0, I.height):
        for col in range(0, I.width):
            red, green, blue = colorsys.hsv_to_rgb(
                                 H_min + H_mul * H[row, col], 
                                 S_min + S_mul * S[row, col], 
                                 V_min + V_mul * V[row, col])
    
            red = int(np.round( red * 255 ))
            green = int(np.round( green * 255 ))
            blue = int(np.round( blue * 255 ))
            P = (red, green, blue)
            I.putpixel((col, row), P)
            
    return I


def mat_to_mapped(A, mp):
    n_rows = A.shape[0]
    n_cols = A.shape[1]
    A, nClrs = flat_index(A)
    mp = resize_color_map(mp, nClrs)*255

    I = IP.new('RGB', (n_cols, n_rows))

    for r in range(0, n_rows):
        for c in range(0, n_cols):
            I.putpixel((c,r), tuple(np.uint8(mp[A[r,c], :])))

    return I





