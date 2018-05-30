
import os
import sys
import time

import numpy as np

from IPython.display import display
import ipywidgets as widgets

sys.path.insert(1, '../src')
import numcolorpy as ncp

def get_primes_glyfic(n, n_start=1, prime_sym=" * ", not_prime_sym='   '):
    """ Usage: glyfic_str = get_primes_glyfic(n, n_start=1, prime_sym=" * ", not_prime_sym='   ') 
    printable string view of pattern of primes
    """
    glyfic_str = ''

    spiral_primes_mask, top_prime, min_prime, num_primes = get_spiral_primes_mask(n, n_start)
    actual_number_of_primes = spiral_primes_mask.sum()
    
    for k in range(0, spiral_primes_mask.shape[0]):
        r = spiral_primes_mask[k,:]
        row_str = ''
        for c in r:
            if c == 1:
                row_str += prime_sym
            else:
                row_str += not_prime_sym

        glyfic_str += row_str + '\n'
        
    return glyfic_str

def show_str_print_primes(n, n_start=1, sym=" * ", not_sym = '   '):
    print(get_primes_glyfic(n, n_start, sym, not_sym))
    
    
def get_prime_number_theorem_estimate(n):
    """ Usage: max_number_of_primes = get_prime_number_theorem_upper_limit(n) """
    
    return np.ceil(n / np.log(n))


def get_spiral_primes_mask(number_rows_cols=7, n_start=1, clockwise=True):
    """ Usage:
    spiral_primes_mask, top_prime, min_prime, num_primes = get_spiral_primes_mask(number_rows_cols, clockwise) 
    """
    spiral_primes_mask = spiral_raster_index(number_rows_cols, n_start, clockwise, just_primes=True)
    top_prime = spiral_primes_mask.max()
    min_prime = spiral_primes_mask.min()
    num_primes = (spiral_primes_mask != 0).sum()
    return spiral_primes_mask.astype(bool).astype(int), top_prime, min_prime, num_primes

def get_prime_divisors_list(x, fname):
    """  """
    max_divisor = x // np.sqrt(x)
    prime_divisors_array = np.loadtxt(fname, dtype=int)
    if x in prime_divisors_array:
        divisors_array = np.array(x)
    else:
        divisors_array = (prime_divisors_array[prime_divisors_array <= max_divisor]).astype(float)
    # option to append list here
    
    return divisors_array

def write_prime_divisors_list(fname, prime_divisors_array):
    """  """
    np.savetxt(fname, prime_divisors_array, fmt='%i')
    
def get_primes_array(max_prime):
    """  """
    array_size = (get_prime_number_theorem_estimate(max_prime) * 1.5).astype(int)
    primes_array = np.zeros(array_size)
    d = 3
    next_loc = 0
    while d < max_prime:
        if is_prime(d):
            primes_array[next_loc] = d
            next_loc += 1
        d += 2
            
    return primes_array[0:next_loc]

def append_primes_array(primes_array, max_prime):
    """  """
    d = primes_array[-1]
    next_loc = max(primes_array.shape)
    if max_prime <= d:
        return primes_array
    array_size = ((get_prime_number_theorem_estimate(max_prime) * 1.5) - next_loc).astype(int)
    appendix_array = np.zeros(array_size)
    apped_array = np.concatenate([primes_array, appendix_array])
    
    while d < max_prime:
        d += 2
        if is_prime(d):
            apped_array[next_loc] = d
            next_loc += 1
            
    return apped_array[0:next_loc-1]

# use the list
# def is_big_prime(x):

def get_next_prime(x, max_search_numbers=1000000):
    """ get the next prime number larger than x """
    if np.mod(x, 2) == 0:
        x += 1
    stopper = 0
    while stopper < max_search_numbers:
        x += 2
        if is_prime(x):
            break
            
    return x


def is_prime(x):
    """ Usage: T_er_F = is_prime(x) """
    if x < 2 or np.mod(x, 2) == 0:
        return False
    sq = x // np.sqrt(x)
    if np.mod(sq, 2) == 0:
        sq -= 1
    for d in range(3, int(sq)+1, 2):
        if np.mod(x, d) == 0:
            return False
    return True


def spiral_raster_index(number_rows_cols=7, n_start=1, clockwise=True, just_primes=False):
    """ Usage: spiral_raster_matrix = spiral_raster_index(number_rows_cols) """
    number_rows_cols_is_even = False
    if np.mod(number_rows_cols, 2) == 0:
        number_rows_cols_is_even = True
        number_rows_cols += 1
        
    spir_ras_mat = np.zeros((number_rows_cols, number_rows_cols))
    r = number_rows_cols // 2
    if clockwise == False:
        rd = 1
    else:
        rd = -1
    c = number_rows_cols // 2
    cd = -1
    rows = [r]
    cols = [c]
    n = max(n_start - 1, 0)
    maxn = n + number_rows_cols**2

    while n < maxn and r >= 0 and r < number_rows_cols:
        for c in cols:
            n += 1
            if just_primes == False or is_prime(n):
                spir_ras_mat[r,c] = n
            
        cd *= (-1)
        c += cd
        cols.append(c)
        cols.reverse()
        
        if c < 0 or c >= number_rows_cols:
            break
        
        for r in rows:
            n += 1
            if just_primes == False or is_prime(n):
                spir_ras_mat[r,c] = n
            
        rd *= (-1)
        r += rd
        rows.append(r)
        rows.reverse()
        
    if number_rows_cols_is_even == True:
        if clockwise == True:
            spir_ras_mat = spir_ras_mat[1:, 1:]
        else:
            spir_ras_mat = spir_ras_mat[:-1, 1:]
            
    return spir_ras_mat.astype(int)    