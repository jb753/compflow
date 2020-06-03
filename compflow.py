"""Perfect gas compressible flow relations.

TODO:
Add per-variable functions to call without branching
"""

import numpy as np
from scipy.optimize import newton

def check_input(y, ga):
    if ga < 1.:
        raise ValueError('Specific heat ratio must be at least unity.')
    if np.any(y < 0.):
        raise ValueError('Input quantity must be positive.')


def to_Ma(var, Y_in, ga, supersonic=False):
    """Invert the Mach number relations by solving iteratively."""

    Y = np.asarray(Y_in)
    check_input(Y, ga)

    # Choose initial guess on sub or supersonic branch
    if supersonic or var in ['Posh_Po', 'Mash']:
        Ma_guess = np.array(1.5)
    else:
        Ma_guess = np.array(0.3)

    Ma_out = np.ones_like(Y)

    # Get indices for non-physical values
    if var == 'mcpTo_APo':
        ich = Y > from_Ma('mcpTo_APo', 1., ga)
    elif var in ['Po_P', 'To_T', 'rhoo_rho', 'A_Acrit']:
        ich = Y < 1.
    elif var in ['Posh_Po', 'Mash']:
        ich = Y > 1.
    else:
        ich = np.full(np.shape(Y), False)

    # Return NaN if the value is not physical
    Ma_out[ich] = np.nan

    # Don't try to solve if all are non-physical
    if np.any(~ich):

        # Check if an explicit inversion exists
        if var == 'To_T':
            Ma_out[~ich] = Ma_from_To_T(Y[~ich], ga)

        elif var == 'Po_P':
            Ma_out[~ich] = Ma_from_Po_P(Y[~ich], ga)

        elif var == 'rhoo_rho':
            Ma_out[~ich] = Ma_from_rhoo_rho(Y[~ich], ga)

        elif var == 'V_cpTo':
            Ma_out[~ich] = Ma_from_V_cpTo(Y[~ich], ga)

        else:

            # Velocity and mass flow functions
            if var == 'V_cpTo':
                f = V_cpTo_from_Ma
                fp = der_V_cpTo_from_Ma

            if var == 'mcpTo_APo':
                f = mcpTo_APo_from_Ma
                fp = der_mcpTo_APo_from_Ma

            if var == 'mcpTo_AP':
                f = mcpTo_AP_from_Ma
                fp = der_mcpTo_AP_from_Ma

            # Choking area
            if var == 'A_Acrit':
                f = A_Acrit_from_Ma
                fp = der_A_Acrit_from_Ma

            # Post-shock Mach
            if var == 'Mash':
                f = Mash_from_Ma
                fp = der_Mash_from_Ma

            # Shock pressure ratio
            if var == 'Posh_Po':
                f = Posh_Po_from_Ma
                fp = der_Posh_Po_from_Ma

            # Wrapper functions for the iterative solve
            def err(Ma_i):
                return f(Ma_i, ga) - Y[~ich]

            def jac(Ma_i):
                return fp(Ma_i, ga)

            # Newton iteration
            Ma_out[~ich] = newton(
                err, np.ones_like(Y[~ich]) * Ma_guess, fprime=jac)

    return Ma_out


def Ma_from_To_T(To_T, ga):
    return np.sqrt((To_T - 1.) * 2. / (ga - 1.))


def Ma_from_Po_P(Po_P, ga):
    return np.sqrt((Po_P ** ((ga - 1.) / ga) - 1.) * 2. / (ga - 1.))


def Ma_from_rhoo_rho(rhoo_rho, ga):
    return np.sqrt((rhoo_rho ** (ga - 1.) - 1.) * 2. / (ga - 1.))


def Ma_from_V_cpTo(V_cpTo, ga):
    return np.sqrt( V_cpTo **2. / (ga - 1.) / (1. - 0.5 * V_cpTo **2.))


def Ma_from_mcpTo_APo(mcpTo_APo, ga, supersonic=False):
    def f(x):
        return mcpTo_APo_from_Ma(x, ga) - mcpTo_APo

    def fp(x):
        return der_mcpTo_APo_from_Ma(x, ga)

    if supersonic:
        Ma_guess = 1.3 * np.ones_like(mcpTo_APo)
    else:
        Ma_guess = 0.3 * np.ones_like(mcpTo_APo)

    return newton(f, Ma_guess, fprime=fp)


def Ma_from_mcpTo_AP(mcpTo_AP, ga):
    def f(x):
        return mcpTo_AP_from_Ma(x, ga) - mcpTo_AP

    def fp(x):
        return der_mcpTo_AP_from_Ma(x, ga)

    Ma_guess = 0.3 * np.ones_like(mcpTo_AP)

    return newton(f, Ma_guess, fprime=fp)


def Ma_from_A_Acrit(A_Acrit, ga, supersonic=False):
    def f(x):
        return A_Acrit_from_Ma(x, ga) - A_Acrit

    def fp(x):
        return der_A_Acrit_from_Ma(x, ga)

    if supersonic:
        Ma_guess = 1.3 * np.ones_like(A_Acrit)
    else:
        Ma_guess = 0.3 * np.ones_like(A_Acrit)

    return newton(f, Ma_guess, fprime=fp)


def Ma_from_Mash(Mash, ga):
    def f(x):
        return Mash_from_Ma(x, ga) - Mash

    def fp(x):
        return der_Mash_from_Ma(x, ga)

    Ma_guess = 1.3 * np.ones_like(Mash)

    return newton(f, Ma_guess, fprime=fp)


def Ma_from_Posh_Po(Posh_Po, ga):
    def f(x):
        return Posh_Po_from_Ma(x, ga) - Posh_Po

    def fp(x):
        return der_Posh_Po_from_Ma(x, ga)

    Ma_guess = 1.3 * np.ones_like(Posh_Po)

    return newton(f, Ma_guess, fprime=fp)


def To_T_from_Ma(Ma, ga):
    return 1. + 0.5 * (ga - 1.0) * Ma ** 2.


def Po_P_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    return To_T ** (ga / (ga - 1.0))


def rhoo_rho_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    return To_T ** (1. / (ga - 1.0))


def V_cpTo_from_Ma(Ma, ga):
    with np.errstate(invalid='ignore'):
        V_cpTo = np.asarray(np.sqrt(
            (ga - 1.0) * Ma * Ma / (1. + (ga - 1.)*Ma*Ma/2.) ))
    V_cpTo[np.isinf(Ma)] = np.sqrt(2.)
    return V_cpTo


def mcpTo_APo_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    with np.errstate(invalid='ignore'):
        mcpTo_APo = np.asarray(ga / np.sqrt(ga - 1.0) * Ma
                            * To_T ** (-0.5 * (ga + 1.0) / (ga - 1.0)))
    mcpTo_APo[np.isinf(Ma)] = 0.
    return mcpTo_APo


def mcpTo_AP_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    return ga / np.sqrt(ga - 1.0) * Ma * To_T ** 0.5


def A_Acrit_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    # Safe reciprocal of Ma
    with np.errstate(divide='ignore'):
        recip_Ma = 1. / Ma
    return (recip_Ma
            * (2. / (ga + 1.0) * To_T) ** (0.5 * (ga + 1.0) / (ga - 1.0)))


def Malimsh_from_ga(ga):
    return np.sqrt((ga - 1.) / ga / 2.) + 0.001


def Mash_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    Mash = np.asarray(np.ones_like(Ma) * np.nan)
    Malimsh = Malimsh_from_ga(ga)
    with np.errstate(invalid='ignore'):
        Mash[Ma >= Malimsh] = (To_T[Ma >= Malimsh]
                / (ga * Ma[Ma >= Malimsh] ** 2. - 0.5 * (ga - 1.0))) ** 0.5
    Mash[np.isinf(Ma)] = np.sqrt( (ga - 1.) / ga / 2.)
    return Mash


def Posh_Po_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    Posh_Po = np.asarray(np.ones_like(Ma) * np.nan)
    Malimsh = Malimsh_from_ga(ga)
    with np.errstate(invalid='ignore'):
        A = 0.5 * (ga + 1.) * Ma ** 2. / To_T
        B = 2. * ga / (ga + 1.0) * Ma ** 2. - 1. / (ga + 1.0) * (ga - 1.0)
        Posh_Po[Ma >= Malimsh] = (A[Ma >= Malimsh] ** (ga / (ga - 1.0))
                                  * B[Ma >= Malimsh] ** (-1. / (ga - 1.)))
    Posh_Po[np.isinf(Ma)] = 0.
    return Posh_Po


def from_Ma(var, Ma_in, ga):
    """Evaluate compressible flow quantities as explicit functions of Ma."""

    Ma = np.asarray(Ma_in)
    check_input(Ma, ga)

    # Simple ratios
    if var == 'To_T':
        return To_T_from_Ma(Ma, ga)

    if var == 'Po_P':
        return Po_P_from_Ma(Ma, ga)

    if var == 'rhoo_rho':
        return rhoo_rho_from_Ma(Ma, ga)

    # Velocity and mass flow functions
    if var == 'V_cpTo':
        return V_cpTo_from_Ma(Ma, ga)

    if var == 'mcpTo_APo':
        return mcpTo_APo_from_Ma(Ma, ga)

    if var == 'mcpTo_AP':
        return mcpTo_AP_from_Ma(Ma, ga)

    # Choking area
    if var == 'A_Acrit':
        return A_Acrit_from_Ma(Ma, ga)

    # Post-shock Mach
    if var == 'Mash':
        return Mash_from_Ma(Ma, ga)

    # Shock pressure ratio
    if var == 'Posh_Po':
        return Posh_Po_from_Ma(Ma, ga)

    # Throw an error if we don't recognise the requested variable
    raise ValueError('Invalid quantity requested: {}.'.format(var))


def der_To_T_from_Ma(Ma, ga):
    return (ga - 1.) * Ma


def der_Po_P_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    return ga * Ma * To_T ** (ga / (ga - 1.) - 1.)


def der_rhoo_rho_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    return Ma * To_T ** (1. / (ga - 1.) - 1.)


def der_V_cpTo_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    return (np.sqrt(ga - 1.)
            * (To_T ** -0.5 - 0.5 * (ga - 1.) * Ma ** 2. * To_T ** -1.5))


def der_mcpTo_APo_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    return ((ga / np.sqrt(ga - 1.)
             * (To_T ** (-0.5 * (ga + 1.) / (ga - 1.))
                 - 0.5 * (ga + 1.) * Ma ** 2.
                 * To_T ** (-0.5 * (ga + 1.) / (ga - 1.) - 1.))))


def der_mcpTo_AP_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    return (ga / np.sqrt(ga - 1.)
            * (To_T ** 0.5 + 0.5 * (ga - 1.) * Ma ** 2. * To_T ** -0.5))


def der_A_Acrit_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    # Safe reciprocal of Ma
    with np.errstate(divide='ignore'):
        recip_Ma = 1. / Ma
    return ((2. / (ga + 1.) * To_T) ** (0.5 * (ga + 1.) / (ga - 1.))
            * (-recip_Ma ** 2. + 0.5 * (ga + 1.) * To_T ** -1.))


def der_Mash_from_Ma(Ma, ga):
    der_Mash = np.asarray(np.ones_like(Ma) * np.nan)
    Malimsh = Malimsh_from_ga(ga)
    To_T = To_T_from_Ma(Ma, ga)
    A = (ga + 1.) ** 2. * Ma / np.sqrt(2.)
    C = ga * (2 * Ma ** 2. - 1.) + 1
    der_Mash[Ma >= Malimsh] = (-A[Ma >= Malimsh] *
                               To_T[Ma >= Malimsh] ** -.5
                               * C[Ma >= Malimsh] ** -1.5)
    return der_Mash


def der_Posh_Po_from_Ma(Ma, ga):
    der_Posh_Po = np.asarray(np.ones_like(Ma) * np.nan)
    Malimsh = Malimsh_from_ga(ga)
    To_T = To_T_from_Ma(Ma, ga)
    A = ga * Ma * (Ma ** 2. - 1.) ** 2. / To_T ** 2.
    B = (ga + 1.) * Ma ** 2. / To_T / 2.
    C = 2. * ga / (ga + 1.) * Ma ** 2. - 1. / (ga + 1.) * (ga - 1.)
    der_Posh_Po[Ma >= Malimsh] = (-A[Ma >= Malimsh]
                                  * B[Ma >= Malimsh] ** (1. / (ga - 1.))
                                  * C[Ma >= Malimsh] ** (-ga / (ga - 1.)))
    return der_Posh_Po


def derivative_from_Ma(var, Ma_in, ga):
    """Evaluate compressible flow quantity derivatives as explict functions """

    Ma = np.asarray(Ma_in)
    check_input(Ma, ga)

    # Simple ratios
    if var == 'To_T':
        return der_To_T_from_Ma(Ma, ga)

    if var == 'Po_P':
        return der_Po_P_from_Ma(Ma, ga)

    if var == 'rhoo_rho':
        return der_rhoo_rho_from_Ma(Ma, ga)

    # Velocity and mass flow functions
    if var == 'V_cpTo':
        return der_V_cpTo_from_Ma(Ma, ga)

    if var == 'mcpTo_APo':
        return der_mcpTo_APo_from_Ma(Ma, ga)

    if var == 'mcpTo_AP':
        return der_mcpTo_AP_from_Ma(Ma, ga)

    # Choking area
    if var == 'A_Acrit':
        return der_A_Acrit_from_Ma(Ma, ga)

    # Post-shock Mack number
    if var == 'Mash':
        return der_Mash_from_Ma(Ma, ga)

    # Shock pressure ratio
    if var == 'Posh_Po':
        return der_Posh_Po_from_Ma(Ma, ga)

    # Throw an error if we don't recognise the requested variable
    raise ValueError('Invalid quantity requested: {}.'.format(var))
