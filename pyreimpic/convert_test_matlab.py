"""
Utilities for conversion of matlab mathematical functions code to python
example:
S_matlab = 'Z = Z^(2*Z^(-c(1)^(Z^-c(2))^(Z^-c(3))^(Z^-c(4))^(Z^-c(5))))'
S_python = matlab_eq_to_py_eq_str(S_matlab, par_char='c', replace_char='p')

returns:
Z = Z**(2*Z**(-p[0]**(Z**-p[1])**(Z**-p[2])**(Z**-p[3])**(Z**-p[4])))

"""


def matlab_eq_to_py_eq_str(S, par_char='c', replace_char='p'):
    """ Usage: S = matlab_eq_to_py_eq_str(S_matlab, par_char='c', replace_char='p')
    Caution: par_char must not be ambuguous in context of whole string

    Args:
        S:              ttring representation of matlab one line equation
        par_char:       parameters variable input name   -   caution
        replace_char:   parameters variable new name

    Returns:
        S_py:           python version of matlab equation with three conversion features:
                            up carret power operator replaced by **
                            parameter variable values reduced by one & parentheses replaced by braces
                            parameter identifier (single) character replaced-renamed

    """
    S_list = list(S)
    for idx in range(len(S_list)):
        if S[idx] == par_char:
            n_start = S[idx:].find('(')
            S_list[idx + n_start] = '['
            n_end = S[idx:].find(')')
            S_list[idx + n_end] = ']'
            S_list[idx + n_start + 1:idx + n_end] = list('%i' % (int(S[idx + n_start + 1:idx + n_end]) - 1))
            S_list[idx] = replace_char

    S_py = ''.join(S_list)

    return S_py.replace('^', '**')