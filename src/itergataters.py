# -*- coding: utf-8 -*-
"""
itergators.py
@author: del
lanier4@illinois.edu
mradmstr514226508@gmail.com
"""
import os
import time
import itertools
import multiprocessing as mp

import numpy as np

import z_plane as zp


def get_primitives(list_tuple, par_set):
    """ ET, Z, Z0 = get_primitives(list_tuple, par_set)
    Args:
        list_tuple:
        par_set:            parameters dictionary with keys: all needed for complex frame,
                            'dir_path', 'it_max', 'max_d'
        delete_temp_dir:    True or False ? delete the temporary directory after assembling Z
    Returns:
        ET:                 Escape Time matrix (may be float - fractional escape times possible)
        Z:                  Complex matrix after iteration
        Z0:                 Complex plane matrix before iteration
    """        
    if 'delete_temp_dir' in par_set:
        delete_temp_dir = par_set['delete_temp_dir']
    else:
        delete_temp_dir = True
    
    par_set['tmp_dir'] = get_tmp_dir(par_set['dir_path'], 'tmp')
    complex_frame, par_set = zp.get_frame_from_dict(par_set)
    n_cores = mp.cpu_count()
    range_enumeration = np.int_(range(0, par_set['n_rows']))
    core_pool = mp.Pool(processes=n_cores)
    core_pool.starmap(write_row, 
                     zip(itertools.repeat(complex_frame),
                         itertools.repeat(list_tuple),
                         itertools.repeat(par_set),
                         range_enumeration))
    core_pool.close()
    core_pool.join()

    Z0, Z, ET = assemble_rows(par_set)
    
    if delete_temp_dir: remove_tmp_dir(par_set['tmp_dir'])

    return ET, Z, Z0


def assemble_rows(par_set):
    """ Z0, Z, ET = assemble_rows(par_set)
        read the temporary files into the output matrices Z0, Z, and ET
    Args:
        par_set:            parameters dictionary with keys:
                            'n_rows', 'n_cols', 'tmp_dir'
    Returns:
        ET:                 Escape Time matrix (may be float - fractional escape times possible)
        Z:                  Complex matrix after iteration
        Z0:                 Complex plane matrix before iteration
    """
    n_rows = par_set['n_rows']
    n_cols = par_set['n_cols']
    tmp_dir = par_set['tmp_dir']
    
    Z0 = np.zeros((n_rows, n_cols),complex)
    Z = np.zeros((n_rows, n_cols),complex)
    ET = np.zeros((n_rows, n_cols))
    
    dir_listing = os.listdir(tmp_dir)
    for f in dir_listing:
        fnm = os.path.join(tmp_dir, f)
        three_sum = np.load(fnm)
        row_number = int(np.imag(three_sum[1, 0]))
        Z0[row_number, :] = three_sum[0, :]
        ET[row_number, :] = np.real(three_sum[1, :])
        Z[row_number, :] = three_sum[2, :]

    return Z0, Z, ET


def write_row(complex_frame, list_tuple, par_set, row_number):
    """ write_row(complex_frame, list_tuple, par_set, row_number)
        Process one row of complex plane in the context of the list_tuple equation definition
        and write the row to the temp directory.

    Args:
        complex_frame:  z_plane complex frame dictionary of points
        list_tuple:     list of tuples data struct:   [(function_handle, ([par1, par2,... parN]) )]
        par_set:        parameters dictionary with keys:
                        'tmp_dir', 'it_max', 'max_d'
        row_number:     integer number of the row to process
    """
    if 'eq_order' in par_set:
        eq_order = par_set['eq_order']
    else:
        eq_order = 0

    if 'RANDOM_PLANE' in par_set and par_set['RANDOM_PLANE'] == True:
        SF = np.abs(complex_frame['upper_left'] - complex_frame['bottom_right'])
        row_array = (np.random.randn(par_set['n_cols']) + np.random.randn(par_set['n_cols']) * 1j) * SF
    else:
        left_style = np.linspace(complex_frame['upper_left'],
                                 complex_frame['bottom_left'], par_set['n_rows'])
        right_style = np.linspace(complex_frame['upper_right'],
                                  complex_frame['bottom_right'], par_set['n_rows'])
        row_array = np.linspace(left_style[row_number], right_style[row_number], par_set['n_cols'])

    it_max = par_set['it_max']
    max_d = par_set['max_d']
    Z_arr = np.zeros((3,row_array.size), complex)

    if eq_order == 0:
        col = 0
        for Z0 in row_array:
            Z_arr[0, col] = Z0
            Z_arr[1, col], Z_arr[2, col] = tuplerator(list_tuple, Z0, it_max, max_d)
            col += 1
    elif eq_order == 2:
        col = 0
        for Z0 in row_array:
            Z_arr[0, col] = Z0
            Z_arr[1, col], Z_arr[2, col] = tuplerator_2(list_tuple, Z0, it_max, max_d)
            col += 1
    elif eq_order == -1:
        col = 0
        for Z0 in row_array:
            Z_arr[0, col] = Z0
            Z_arr[1, col], Z_arr[2, col] = tuplerator_3(list_tuple, Z0, it_max, max_d, par_set)
            col += 1

    Z_arr[1, :] = Z_arr[1, :] + complex(0.0, float(row_number)) # row number as imaginary part
    file_name = os.path.join(par_set['tmp_dir'], ahora_seq_name('row_%d_'%(row_number), '.txt'))
    with open(file_name, 'wb') as file_handle:
        Z_arr.dump(file_handle)

        
def tuplerator_3(list_tuple, Z0, it_max, max_d, par_set):
    """ ET, Z = tuplerator(list_tuple, Z0, it_max, max_d)
        function iterator for two tuple list (with at least one tuple)
        list_tuple = [(function_1, (p1, p2,... , pn)), (function_2, (p1, p2,... , pn))]
        
    Args:
        list_tuple: [ ( function_handle_01,(fh_01_p1, fh_01_p2,...) ),... (fh_n,(p1, p2,...)) ]
                    list of tuples: function handle (not lambda fcn) & internal parameters
        Z0:         complex vector - starting point
        it_max:     maximum number of iterations
        max_d:      quit iterating and return distance
        
    Returns:
        ET:         number of iterations to finishing point
        Z:          complex vector - finishing point
        
    Prototype:
        def jsV(Z, P):
            Z = Z**P[0] - P[1]
            return Z
    """
    ET = 0
    Z = Z0
    d = 0
    while (ET <= it_max) & (np.isfinite(d)) & (d < max_d):
        ET += 1
        Z_was = Z
        try:
            for fcn_hndl, P in list_tuple:
                Z = fcn_hndl(Z, P, Z0, ET, d, par_set)
                # checkered_daemon(Z, p, Z0=None, ET=None, d=0, par_set=None)
            d = np.abs(Z - Z0)
        except:
            return max(1, ET - 1), Z_was
    # if (ET >= it_max) & (np.isfinite(d)) & (d < max_d):    
    if (np.isfinite(d)):
        return ET, Z
    else:
        return max(1, ET - 1), Z_was
        
def tuplerator_2(list_tuple, Z0, it_max, max_d):
    ET = 0
    Z = Z0
    Zm1 = Z0
    Zm2 = Z0
    d = 0
    while (ET <= it_max) & (np.isfinite(d)) & (d < max_d):
        Z_was = Z
        try:
            for fcn_hndl, P in list_tuple:
                Z, Zm1, Zm2 = fcn_hndl(Z, P, Z0, ET, Zm1, Zm2)
            d = np.abs(Z - Z0)

        except:
            return max(1, ET - 1), Z_was
        ET += 1
    # if (ET >= it_max) & (np.isfinite(d)) & (d < max_d):
    if (np.isfinite(d)):
        return ET, Z
    else:
        return max(1, ET - 1), Z_was

def tuplerator(list_tuple, Z0, it_max, max_d):
    """ ET, Z = tuplerator(list_tuple, Z0, it_max, max_d)
        function iterator for two tuple list (with at least one tuple)
        list_tuple = [(function_1, (p1, p2,... , pn)), (function_2, (p1, p2,... , pn))]
        
    Args:
        list_tuple: [ ( function_handle_01,(fh_01_p1, fh_01_p2,...) ),... (fh_n,(p1, p2,...)) ]
                    list of tuples: function handle (not lambda fcn) & internal parameters
        Z0:         complex vector - starting point
        it_max:     maximum number of iterations
        max_d:      quit iterating and return distance
        
    Returns:
        ET:         number of iterations to finishing point
        Z:          complex vector - finishing point
        
    Prototype:
        def jsV(Z, P):
            Z = Z**P[0] - P[1]
            return Z
    """
    ET = 0
    Z = Z0
    d = 0
    while (ET <= it_max) & (np.isfinite(d)) & (d < max_d):
        ET += 1
        Z_was = Z
        try:
            for fcn_hndl, P in list_tuple:
                Z = fcn_hndl(Z, P, Z0)
            d = np.abs(Z - Z0)
        except:
            return max(1, ET - 1), Z_was
    # if (ET >= it_max) & (np.isfinite(d)) & (d < max_d):    
    if (np.isfinite(d)):
        return ET, Z
    else:
        return max(1, ET - 1), Z_was
    
    
    
    
def ahora_seq_name(prefi_str=None, suffi_str=None):
    """ alpha_time_stamped_name = ahora_seq_name(prefi_str, suffi_str)
        locally unique (1/1e12 sec) alphanumeric-integer time stamp name
    
    Args:
        n_digits:  number of digits accuracy in the decimal seconds part of the name 
        prefi_str: prepend to time_stamp
        suffi_str: postpent to time stamp
        
    Returns:
        now_name:  a time stamped file name with prefix and suffix as input 
    """
    time_decimal_shift = 1e12
    ahora_nombre = '%d'%(time.time() * time_decimal_shift)
    if prefi_str is None: prefi_str = ''
    if suffi_str is None: suffi_str = ''

    return prefi_str + ahora_nombre + suffi_str


def int_to_alpha(N):
    """ alpha_N = int_to_alpha(N)
        convert integer to spreadsheet column lettering

    Args:    an integer representation of the letters like 136964

    Returns: a string of capital letters like 'GTOV'
    """
    B = 26
    ALPH = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    max_mag = 100
    M = 0

    while (N > B ** M) & (M <= max_mag):
        M += 1

    S = []
    for nada in range(0, M):
        M -= 1
        D = np.int_(np.floor(N / B ** M))
        N = N - D * B ** M
        S.append(ALPH[D - 1])

    alpha_N = ''.join(S)

    return alpha_N


def seq_time_id_string(time_units=1e6):
    """ alpha_N = seq_time_id_string(time_units)
        current time in time_units as an integer encoded in a base 26 caps string

    Args:
        time_units:  1e6 is micro seconds, 1 <= time_units <= 1e9

    Returns:
        tid_str:     current time in time_units as an integer encoded in a base 26 caps string
    """
    return int_to_alpha(int(time.time() * np.maximum(np.minimum(time_units, 1e9), 1)))


def get_tmp_dir(dir_path, dir_name, timestamp=None):
    """ new_dir_name = get_tmp_dir(dir_path, dir_name, timestamp)
        create a "dir_name" with time stamp directory

    Args:
        dir_path:       an existing directory such as the run directory.
        dir_name:       new directory name to add to dir_path directory.
        timestamp:      optional - if not input a microsecond stamp will be added.
    Returns:
        new_dir_name:   name of directory just created
    """
    if timestamp is None:
        timestamp = seq_time_id_string()
    new_dir_name = os.path.join(dir_path, dir_name + timestamp)
    os.mkdir(new_dir_name)

    return new_dir_name


def remove_tmp_dir(dir_name):
    """ remove_tmp_dir(dir_name)
        remove the directory and all the files it contains.

    Args:
        dir_name: name of a directory with no sub-directories.
    """
    try:
        dir_list = os.listdir(dir_name)
        if len(dir_list) > 0:
            for file_name in dir_list:
                os.remove(os.path.join(dir_name, file_name))
        os.rmdir(dir_name)
    except:
        print('function remove_tmp_dir:\n%s\n(directory DNE)' % (dir_name))
        pass

    return


class ETA_iterator:

    def __init__(self, it_max=32, max_d=12):
        self._it_max = it_max
        self._max_d = max_d
