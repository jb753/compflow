"""Functions to calculate one-dimensional compressible flow quantities.

TODO:
Add per-variable functions to call without branching
"""

import numpy as np
from scipy.optimize import newton


def to_Ma(var, Y_in, ga, supersonic=False, validate=True):
    """Invert the Mach number relations by solving iteratively."""

    if validate:
        # Validate input data
        Y = np.asarray_chkfinite(Y_in)
        if ga <= 0.:
            raise ValueError('Specific heat ratio must be positive.')
        if np.any(Y < 0.):
            raise ValueError('Input quantity must be positive.')
    else:
        Y = np.asarray(Y_in)

    # Choose initial guess on sub or supersonic branch
    if supersonic or var in ['Posh_Po', 'Mash']:
        Ma_guess = np.array(1.5)
    else:
        Ma_guess = np.array(0.3)

    Ma_out = np.ones_like(Y)

    # Get indices for non-physical values
    if var == 'mcpTo_APo':
        ich = Y > from_Ma('mcpTo_APo', 1., ga, validate=False)
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
            Ma_out[~ich] = np.sqrt((Y[~ich] - 1.) * 2. / (ga - 1.))

        elif var == 'Po_P':
            Ma_out[~ich] = np.sqrt((Y[~ich] ** ((ga - 1.) / ga) - 1.) * 2. / (ga - 1.))

        elif var == 'rhoo_rho':
            Ma_out[~ich] = np.sqrt((Y[~ich] ** (ga - 1.) - 1.) * 2. / (ga - 1.))

        else:

            # Wrapper functions for the iterative solve
            def err(Ma_i):
                return from_Ma(var, Ma_i, ga, validate=False) - Y[~ich]

            def jac(Ma_i):
                return derivative_from_Ma(var, Ma_i, ga, validate=False)

            # Newton iteration
            Ma_out[~ich] = newton(err, np.ones_like(Y[~ich]) * Ma_guess, fprime=jac)

    return Ma_out


def from_Ma(var, Ma_in, ga_in, validate=True):
    """Evaluate compressible flow quantities as explicit functions of Ma."""

    if validate:
        ga = np.asarray_chkfinite(ga_in)
        Ma = np.asarray_chkfinite(Ma_in)
        if ga <= 0.:
            raise ValueError('Specific heat ratio must be positive.')
        if np.any(Ma < 0.):
            raise ValueError('Mach number must be positive.')
    else:
        ga = np.asarray(ga_in)
        Ma = np.asarray(Ma_in)

    # Stagnation temperature ratio appears in every expression
    To_T = 1. + 0.5 * (ga - 1.0) * Ma ** 2.

    # Simple ratios
    if var == 'To_T':
        return To_T

    if var == 'Po_P':
        return To_T ** (ga / (ga - 1.0))

    if var == 'rhoo_rho':
        return To_T ** (1. / (ga - 1.0))

    # Velocity and mass flow functions
    if var == 'V_cpTo':
        return np.sqrt(ga - 1.0) * Ma * To_T ** -0.5

    if var == 'mcpTo_APo':
        return ga / np.sqrt(ga - 1.0) * Ma * To_T ** (-0.5 * (ga + 1.0) / (ga - 1.0))

    if var == 'mcpTo_AP':
        return ga / np.sqrt(ga - 1.0) * Ma * To_T ** 0.5

    # Choking area
    if var == 'A_Acrit':
        # Safe reciprocal of Ma
        with np.errstate(divide='ignore'):
            recip_Ma = 1. / Ma
        return recip_Ma * (2. / (ga + 1.0) * To_T) ** (0.5 * (ga + 1.0) / (ga - 1.0))

    Malimsh = np.sqrt((ga - 1.) / ga / 2.) + 0.001  # Mathematical limit on pre-shock Ma

    # Post-shock Mach
    if var == 'Mash':
        Mash = np.asarray(np.ones_like(Ma) * np.nan)
        Mash[Ma >= Malimsh] = (To_T[Ma >= Malimsh] / (ga * Ma[Ma >= Malimsh] ** 2. - 0.5 * (ga - 1.0))) ** 0.5
        return Mash

    # Shock pressure ratio
    if var == 'Posh_Po':
        Posh_Po = np.asarray(np.ones_like(Ma) * np.nan)
        A = 0.5 * (ga + 1.) * Ma ** 2. / To_T
        B = 2. * ga / (ga + 1.0) * Ma ** 2. - 1. / (ga + 1.0) * (ga - 1.0)
        Posh_Po[Ma >= Malimsh] = A[Ma >= Malimsh] ** (ga / (ga - 1.0)) * B[Ma >= Malimsh] ** (-1. / (ga - 1.))
        return Posh_Po

    # Throw an error if we don't recognise the requested variable
    raise ValueError('Invalid quantity requested: {}.'.format(var))


def derivative_from_Ma(var, Ma_in, ga_in, validate=True):
    """Evaluate compressible flow quantity derivatives as explict functions of Ma."""

    if validate:
        ga = np.asarray_chkfinite(ga_in)
        Ma = np.asarray_chkfinite(Ma_in)
        if ga <= 0.:
            raise ValueError('Specific heat ratio must be positive.')
        if np.any(Ma < 0.):
            raise ValueError('Mach number must be positive.')
    else:
        ga = np.asarray(ga_in)
        Ma = np.asarray(Ma_in)

    # Stagnation temperature ratio appears in every expression
    To_T = 1. + 0.5 * (ga - 1.) * Ma ** 2.

    # Simple ratios
    if var == 'To_T':
        return (ga - 1.) * Ma

    if var == 'Po_P':
        return ga * Ma * To_T ** (ga / (ga - 1.) - 1.)

    if var == 'rhoo_rho':
        return Ma * To_T ** (1. / (ga - 1.) - 1.)

    # Velocity and mass flow functions
    if var == 'V_cpTo':
        return np.sqrt(ga - 1.) * (To_T ** -0.5 - 0.5 * (ga - 1.) * Ma ** 2. * To_T ** -1.5)

    if var == 'mcpTo_APo':
        return ga / np.sqrt(ga - 1.) * \
               (To_T ** (-0.5 * (ga + 1.) / (ga - 1.))
                - 0.5 * (ga + 1.) * Ma ** 2. * To_T ** (-0.5 * (ga + 1.) / (ga - 1.) - 1.))

    if var == 'mcpTo_AP':
        return ga / np.sqrt(ga - 1.) * (To_T ** 0.5 + 0.5 * (ga - 1.) * Ma ** 2. * To_T ** -0.5)

    # Choking area
    if var == 'A_Acrit':
        # Safe reciprocal of Ma
        with np.errstate(divide='ignore'):
            recip_Ma = 1. / Ma
        return (2. / (ga + 1.) * To_T) ** (0.5 * (ga + 1.) / (ga - 1.)) * (
                -recip_Ma ** 2. + 0.5 * (ga + 1.) * To_T ** -1.)

    # Define shorthand gamma combinations
    Malimsh = np.sqrt(0.5 * (ga - 1.) / ga) + 0.001  # Limit when denominator goes negative

    # Post-shock Mack number
    if var == 'Mash':
        der_Mash = np.asarray(np.ones_like(Ma) * np.nan)
        A = (ga + 1.) ** 2. * Ma / np.sqrt(2.)
        C = ga * (2 * Ma ** 2. - 1.) + 1
        der_Mash[Ma >= Malimsh] = -A[Ma >= Malimsh] * To_T[Ma >= Malimsh] ** -.5 * C[Ma >= Malimsh] ** -1.5
        return der_Mash

    # Shock pressure ratio
    if var == 'Posh_Po':
        der_Posh_Po = np.asarray(np.ones_like(Ma) * np.nan)
        A = ga * Ma * (Ma ** 2. - 1.) ** 2. / To_T ** 2.
        B = (ga + 1.) * Ma ** 2. / To_T / 2.
        C = 2. * ga / (ga + 1.) * Ma ** 2. - 1. / (ga + 1.) * (ga - 1.)
        der_Posh_Po[Ma >= Malimsh] = -A[Ma >= Malimsh] * B[Ma >= Malimsh] ** (1. / gm1) * C[Ma >= Malimsh] ** (
                -ga / (ga - 1.))
        return der_Posh_Po

    # Throw an error if we don't recognise the requested variable
    raise ValueError('Invalid quantity requested: {}.'.format(var))
