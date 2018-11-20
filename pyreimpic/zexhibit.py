# -*- coding: utf-8 -*-
"""
lanier4@illinois.edu
mradmstr514226508@gmail.com

functionality to display & print complex plane
"""
import numpy as np
import matplotlib.pyplot as plt

# print_order controls display of frame in complex_frame_dict_to_string
print_order = ['upper_left', 'top_center', 'upper_right',
               'left_center', 'center_point', 'right_center',
              'bottom_left', 'bottom_center', 'bottom_right']

def complex_to_string(z, N_DEC=6):
    """ format single complex number to string with n decimal places """
    MAX_DEC = 17
    MIN_DEC = 1

    n = max(min(MAX_DEC, round(N_DEC)), MIN_DEC)
    fs = '%%0.%df'%n

    zr = np.real(z)
    zi = np.imag(z)

    if np.sign(zi) < 0:
        s1 = ' ' + fs % zi + 'j'
    else:
        s1 = ' +' + fs % zi + 'j'

    if np.sign(zr) < 0:
        z_str = fs % zr + s1
    else:
        z_str = ' ' + fs % zr + s1

    return z_str


def show_complex_matrix(Z0, N_DEC=3):
    """ display a complex matrix or array """
    SPC = ' ' * 2
    if Z0.shape[0] == Z0.size:
        row_str = ''
        for col in range(0, Z0.shape[0]):
            row_str += complex_to_string(Z0[col], N_DEC) + SPC + '\n'
        print(row_str)

    else:
        for row in range(0, Z0.shape[0]):
            row_str = ''
            for col in range(0, Z0.shape[1]):
                row_str += complex_to_string(Z0[row, col], N_DEC) + SPC

            print(row_str)


def get_aligned_dict_string(d, N_DEC=3):
    """ pretty_string = z_plane.get_aligned_dict_string(d, N_DEC=3) """
    INDENT = 16
    out_string = ''
    for k in sorted(list(d.keys())):

        v = d[k]
        if type(v) == str:
            s = v
        elif v == 0:
            s = '0'
        elif np.iscomplex(v):
            s = complex_to_string(v, N_DEC)
        elif np.round(v) == v:
            s = '%d' % (v)
        else:
            f_str = '%s%s%d%s' % ('%', '0.', N_DEC, 'f')
            s = f_str % (v)

        if len(out_string) == 0:
            out_string = ' ' * max(0, (INDENT - len(k))) + k + ': ' + s
        else:
            out_string = out_string + '\n' + ' ' * max(0, (INDENT - len(k))) + k + ': ' + s

    return out_string + '\n'


def complex_frame_dict_to_string(frame_dict, N_DEC=4):
    """ get a formatted list of strings """
    STR_L = 14
    frame_string = ''
    row = 0
    for k in print_order:
        z_str = complex_to_string(frame_dict[k], N_DEC)
        PAD = ' ' * (STR_L - len(z_str))
        frame_string += k + ':' + PAD + z_str
        row += 1
        if np.mod(row, 3) == 0:
            frame_string += '\n'
        else:
            frame_string += '\t'

    return frame_string


# plot_Z1_Z2()
# plot_polar()


def pyplot_raw_complex(Z, title='Z complex', ylabel='Imaginary', xlabel='Real'):
    """ fig, ax = pyplot_raw_complex(Z, title, ylabel, xlabel)
    plot complex vectors from origin """
    fig = plt.figure()
    ax = plt.axes()

    if isinstance(Z, np.ndarray) or isinstance(Z, list):
        Z_arr = mat2arr(Z)
        for z_ix in range(Z_arr.shape[0]):
            plt.plot([0, Z_arr[z_ix].real],
                     [0, Z_arr[z_ix].imag],
                     'go-', label='python')
        limit = np.max(np.ceil(np.absolute(Z_arr)))

    else:
        plt.plot([0, Z.real], [0, Z.imag], 'bo-', label='python')
        limit = np.max(np.ceil(np.absolute(Z)))

    plt.xlim((-limit, limit))
    plt.ylim((-limit, limit))
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)

    return fig, ax


def mat2arr(M):
    """ unroll a tensor to an array for sequential plotting etc. -- copy in to out """
    if len(M.shape) > 1:
        A = (M + 0.0).ravel()
    else:
        A = M + 0.0
    return A


def pyplot_on_grid_complex(Z0, Z1, title='Z complex', ylabel='height', xlabel='width'):
    """ fig, ax = pyplot_grid_complex(Z0, Z1,title, ylabel, xlabel) """
    Z_start = mat2arr(Z0)
    Z_end = mat2arr(Z1)
    if not isinstance(Z_start, np.ndarray) or isinstance(Z_start, list):
        fig, ax = pyplot_raw_complex(Z, title, ylabel, xlabel)
    else:
        fig = plt.figure()
        ax = plt.axes()
        for arr_idx in range(Z_start.shape[0]):
            plt.plot([Z_start[arr_idx].real, Z_end[arr_idx].real],
                     [Z_start[arr_idx].imag, Z_end[arr_idx].imag],
                     'bo-')

    limit = max(np.max(np.ceil(np.absolute(Z_start))), np.max(np.ceil(np.absolute(Z_end))))

    plt.xlim((-limit, limit))
    plt.ylim((-limit, limit))
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)

    return fig, ax

