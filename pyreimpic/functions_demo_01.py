# -*- coding: utf-8 -*-
"""
lanier4@illinois.edu
mradmstr514226508@gmail.com

first collection of functions
Note that third parameter (ET) is ghosted here because some other functions wannit
"""
import numpy as np


def starfish_ish(Z, p=None, Z0=None, ET=None):
    """
    par_set['zoom'] = 5/8

    Args:
        Z:    a real or complex number
        p:    a real of complex number
    Returns:
        Z:    the result (complex)
    """
    if p is None:
        p = -0.040431211565 + 0.388620268274j
    return Z**(-np.exp(Z**p)**(np.exp(Z**p)**(-np.exp(Z**p)**(np.exp(Z**p)**(-np.exp(Z**p))))))


def starfish_ish_II(Z, p=None, Z0=None, ET=None):
    """
    par_set['zoom'] = 5/8

    Args:
        Z:    a real or complex number
        p:    a real of complex number
    Returns:
        Z:    the result (complex)
    """
    if p is None:
        p = -0.051448293230 + 0.304348945637j
    Z = Z**(-np.exp(Z**p)**(np.exp(Z**p)**(-np.exp(Z**p)**(np.exp(Z**p)**(\
                    -np.exp(Z**p)**(np.exp(Z**p)**(-np.exp(Z**p))))))))
    return Z


def decPwrAFx(Z, p=None, Z0=None, ET=None):
    """
    par_set['zoom'] = 1/8

    Args:
        Z:    a real or complex number
        p:    array real of complex number
    Returns:
        Z:    the result (complex)
    Z = 1/Z - Z^(n*Z^(P(n)^n) / sqrt(pi));
    """
    if p is None:
        p = [np.sqrt(np.pi), 1.13761386, -0.11556857]
    for n in range(1,len(p)):
        Z = 1/Z - Z**(n * Z**(p[n]**n) / p[0])
    return Z


def dreadSkull(Z, p=None, Z0=None, ET=None):
    """
    par_set['zoom'] = 0.4

    Args:
        Z:    a real or complex number
        p:    a real of complex number
    Returns:
        Z:    the result (complex)
        p[0]

    MATLAB:
    Z = (-Z)^(-exp(Z^x)^(exp(Z^x)^(-exp(Z^x)^(exp(Z^x)^(-exp(Z^x)^(exp(Z^x)^(-exp(Z^x))))))))
    """
    if p is None:
        p = -0.295887110004
    ZEP = np.exp(Z ** p)
    Zout = (-Z) ** (-ZEP ** (ZEP ** (-ZEP ** (ZEP ** (-ZEP ** (ZEP ** (-ZEP)))))))
    return Zout

def IslaLace(Z, p=None, Z0=None, ET=None):
    """
    par_set['zoom'] = 1/4

    Args:
        Z:    a real or complex number
        p:    a real of complex number
    Returns:
        Z:    the result (complex)
    """
    if p is None:
        p = [0.444476893762, 0.508164683992 + 0.420921535772j]
    x = p[0]
    c = p[1]
    Z = ( Z**(-x**(Z**(-c))) + x**(-Z**(-c**Z))) * (c**(-Z**(-x**Z)) - Z**(-x**(-Z**(-c))) ) + \
        ( Z**(-x**(Z**(-c))) - x**(-Z**(-c**Z))) * (c**(-Z**(-x**Z)) + Z**(-x**(-Z**(-c))) )
    return Z

def RoyalZ(Z, p=None, Z0=None, ET=None):
    """
    par_set['zoom'] = 1/3

    Args:
        Z:    a real or complex number
        p:    a real of complex number
    Returns:
        Z:    the result (complex)
    """
    if p is None:
        p = [0.340859990388, 0.269282250320, -0.255017720861]
    nc = len(p)
    for n in range(0, nc):
        Z = Z**(-1*np.exp(Z*p[n]))
    return Z

def T_Spake_Z(Z, p, Z0=None, ET=None):
    """ Z = T_Spake_Z(Z, p)
    par_set['zoom'] = 1/3

    Args:
        Z:    a real or complex number
        p:    a real of complex number

    Returns:
        Z:    the result (complex)
    """
    if p is None:
        p = [1.92846051108342, 2.27919841968635, 3.37327534248407, 2.17984103218476]
    d = np.abs(Z-Z0)
    Zxy = np.sqrt(Z/np.abs(Z))
    x = np.real(Zxy)
    y = np.imag(Zxy)*1j
    Z = Z - ( p[0]*x**3 + 3*p[1]*x**2*y + 3*p[2]*x*y**2 + p[3]*y**3 )**(Z*d)

    return Z


def ItchicuPpwrF(Z, p=None, Z0=None, ET=None, Zm1=0, Zm2=0):
    """
    par_set['zoom'] = 0.16

    Args:
        Z:    a real or complex number
        p:    a real of complex number
    Returns:
        Z:    the result (complex)
    """
    if p is None:
        p = [0.56890021, -0.25564542, -0.37746896, -0.29588711, -1.47513451, -0.23400405, 0.11844484]
    for n in range(0, len(p) - 1):
        try:
            Zn = Z ** (2 * Z ** (-(p[n]) ** (Z ** (-p[n + 1]))))
        except:
            return Z
            pass

        if np.isfinite(Zn):
            Z = Zn
        else:
            return Z

    return Z