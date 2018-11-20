# -*- coding: utf-8 -*-
"""
lanier4@illinois.edu
mradmstr514226508@gmail.com

first collection of functions
Note that third parameter (ET) is ghosted here because some other functions wannit
"""
import numpy as np


def starfish_ish(Z, p, Z0=None, ET=None):
    """ Z = starfish_ish(Z, p, Z0=None, ET=None)
    Known parameter:    p = -0.040431211565+0.388620268274j
    Args:
        Z:    a real or complex number
        p:    a real of complex number
    Returns:
        Z:    the result (complex)
    """
    return Z**(-np.exp(Z**p)**(np.exp(Z**p)**(-np.exp(Z**p)**(np.exp(Z**p)**(-np.exp(Z**p))))))


def starfish_ish_II(Z, p, Z0=None, ET=None):
    """ Z = starfish_ish_II(Z, p, Z0=None, ET=None)
    Known parameter:    p = -0.051448293230 + 0.304348945637i
    Args:
        Z:    a real or complex number
        p:    a real of complex number
    Returns:
        Z:    the result (complex)
    """
    Z = Z**(-np.exp(Z**p)**(np.exp(Z**p)**(-np.exp(Z**p)**(np.exp(Z**p)**(\
                    -np.exp(Z**p)**(np.exp(Z**p)**(-np.exp(Z**p))))))))
    return Z


def decPwrAFx(Z, p, Z0=None, ET=None):
    """ Z = decPwrAFx(Z, p, Z0=None, ET=None)
    Known parameter:    p = [np.sqrt(np.pi), 1.13761386, -0.11556857]
    Args:
        Z:    a real or complex number
        p:    array real of complex number
    Returns:
        Z:    the result (complex)
    Z = 1/Z - Z^(n*Z^(P(n)^n) / sqrt(pi));
    """
    for n in range(1,len(P)):
        Z = 1/Z - Z**(n * Z**(P[n]**n) / p[0]);
    return Z


def dreadSkull(Z, p, Z0=None, ET=None):
    """ Z = dreadSkull(Z, p, Z0=None, ET=None)
    Known parameter:    p = -0.295887110004
    Args:
        Z:    a real or complex number
        p:    a real of complex number
    Returns:
        Z:    the result (complex)
        p[0]

    MATLAB:
    Z = (-Z)^(-exp(Z^x)^(exp(Z^x)^(-exp(Z^x)^(exp(Z^x)^(-exp(Z^x)^(exp(Z^x)^(-exp(Z^x))))))))
    """
    ZEP = np.exp(Z ** p)
    Zout = (-Z) ** (-ZEP ** (ZEP ** (-ZEP ** (ZEP ** (-ZEP ** (ZEP ** (-ZEP)))))))
    return Zout

