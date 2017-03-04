"""
complex plane functions and class for getting pixels in three different measuring systems
"""
import numpy as np

print_order = ['upper_left', 'top_center', 'upper_right',
               'left_center', 'center_point', 'right_center',
              'bottom_left', 'bottom_center', 'bottom_right']


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


def get_complex_frame_string(frame_dict, N_DEC=4):
    """ get a formatted list of strings """

    Z_INDENT = 14
    CELL_SPC = 3

    frame_string = []
    row = 0
    for k in print_order:
        z_str = complex_to_string(frame_dict[k], N_DEC)
        spc1 = ' ' * (Z_INDENT - len(k))
        row += 1
        if np.mod(row, 3) == 1:
            row_string = [k + spc1 + z_str]
        elif np.mod(row, 3) == 0:
            spc0 = ' ' * CELL_SPC
            row_string.append(spc0 + k + spc1 + z_str)
            frame_string.append(row_string)
        else:
            spc0 = ' ' * CELL_SPC
            row_string.append(spc0 + k + spc1 + z_str)

    return frame_string


def show_complex_frame(frame_dict):
    """ display a frame_dict in 3 x 3 table format """
    frame_string = get_complex_frame_string(frame_dict)
    for l in frame_string:
        print(l[0], l[1], l[2])



















































""" wolf wolf """