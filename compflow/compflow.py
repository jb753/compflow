"""Perfect gas compressible flow relations.

This module contains functions to convert back and forth between Mach number
and various other non-dimensional flow quantities.

Two interfaces are available. `to_Ma` and `from_Ma` take a string argument that
selects the flow quantity to use, and perform some sanity checks on input data.
There are also lower-level functions which assume the input data is a sensible
numpy array and require no branching, the idea being that these are slightly
faster.

JB June 2020
"""

import numpy as np
from scipy.optimize import newton
from scipy.interpolate import UnivariateSpline

# Initialise empty module-level cache for lookup tables
cache = {}


def generate_lookup(var, ga, atol=1e-7):
    """Generate a lookup table for faster inversions to Mach number.

    Args:
        var (str): Name of flow quantity to be used as independent variable.
        ga (float): Ratio of specific heats.
        atol (float): Absolute tolerance on Mach number, default 1e-7.

    Returns:
        f (InterpolatedUnivariateSpline): Callable interpolator that returns
            Mach number as a function of the requested flow quantity."""

    # Start with a small table, uniformly sampled
    Nk = 100
    y = np.linspace(0., 1., Nk)
    x = from_Ma(var, y, ga)

    # Loop until we reach tolerance
    err_max = np.inf
    N_max = 50
    for n in range(N_max):

        # Make interpolator
        f = UnivariateSpline(x,y,k=1,s=0.,ext='raise')

        # Compute error in Mach at midpoints
        ym = 0.5 * (y[:-1] + y[1:])
        xm = from_Ma(var, ym, ga)
        err = np.abs(f(xm) - ym)
        err_max = np.max(err)

        # Break out of loop if this is good enough
        if err_max < atol:
            break

        # Find indices where err exceeds tolerance and add to table
        x = np.sort(np.hstack((x,xm[err > atol])))
        y = np.sort(np.hstack((y,ym[err > atol])))

    return f


def check_input(y, ga):
    if ga < 1.:
        raise ValueError('Specific heat ratio must be at least unity.')
    if np.any(y < 0.):
        raise ValueError('Input quantity must be positive.')

def get_invalid(var, Y, ga):
    """Return indices for non-physical values."""

    if var == 'mcpTo_APo':
        ich = Y > mcpTo_APo_from_Ma(1., ga)
    elif var in ['Po_P', 'To_T', 'rhoo_rho', 'A_Acrit']:
        ich = Y < 1.
    elif var in ['Posh_Po', 'Mash']:
        ich = Y > 1.
    elif var in ['V_cpTo']:
        ich = Y > np.sqrt(2)
    elif var in ['Mash']:
        ich = Y < np.sqrt((ga - 1.)/ 2. / ga)
    else:
        ich = np.full(np.shape(Y), False)

    return ich

def to_Ma(var, Y_in, ga, supersonic=False, use_lookup=True):
    """Invert the Mach number relations by solving iteratively."""

    # Check if a lookup table exists
    if use_lookup and not supersonic \
        and not var in ['To_T', 'Po_P', 'rhoo_rho', 'V_cpTo']:
        if ga not in cache:
            cache[ga] = {}
        if var not in cache[ga]:
            cache[ga][var] = generate_lookup(var, ga)
        try:
            Ma = cache[ga][var](Y_in)
            return Ma

        except ValueError:
            pass

    # Coerce input to at least a 1D numpy array, so we can use logical indexing
    # and maths operations on it without thinking
    Y = np.atleast_1d(Y_in)

    check_input(Y, ga)

    # nan indicates invalid or non-physical input
    Ma_out = np.empty_like(Y) * np.nan

    # Don't try to solve if all are non-physical
    iv = ~get_invalid(var, Y, ga)
    if np.any(iv):

        # Check if an explicit inversion exists
        if var == 'To_T':
            Ma_out[iv] = Ma_from_To_T(Y[iv], ga)

        elif var == 'Po_P':
            Ma_out[iv] = Ma_from_Po_P(Y[iv], ga)

        elif var == 'rhoo_rho':
            Ma_out[iv] = Ma_from_rhoo_rho(Y[iv], ga)

        elif var == 'V_cpTo':
            Ma_out[iv] = Ma_from_V_cpTo(Y[iv], ga)

        elif var == 'mcpTo_APo':
            Ma_out[iv] = Ma_from_mcpTo_APo(Y[iv], ga, supersonic)

        elif var == 'mcpTo_AP':
            Ma_out[iv] = Ma_from_mcpTo_AP(Y[iv], ga)

        elif var == 'A_Acrit':
            Ma_out[iv] = Ma_from_A_Acrit(Y[iv], ga, supersonic)

        elif var == 'Mash':
            Ma_out[iv] = Ma_from_Mash(Y[iv], ga)

        # Shock pressure ratio
        elif var == 'Posh_Po':
            Ma_out[iv] = Ma_from_Posh_Po(Y[iv], ga)

        else:
            raise ValueError('Bad flow quantity requested.')


    # Return to a scalar if required
    if np.size(Ma_out)==1:
        Ma_out = float(Ma_out)

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
        Ma_guess = 1.5* np.ones_like(mcpTo_APo)
    else:
        Ma_guess = 0.5* np.ones_like(mcpTo_APo)
    return newton(f, Ma_guess , fprime=fp)


def Ma_from_mcpTo_AP(mcpTo_AP, ga):
    def f(x):
        return mcpTo_AP_from_Ma(x, ga) - mcpTo_AP
    def fp(x):
        return der_mcpTo_AP_from_Ma(x, ga)
    Ma_guess = 0.5 * np.ones_like(mcpTo_AP)
    return newton(f, Ma_guess, fprime=fp)


def Ma_from_A_Acrit(A_Acrit, ga, supersonic=False):
    def f(x):
        return A_Acrit_from_Ma(x, ga) - A_Acrit
    def fp(x):
        return der_A_Acrit_from_Ma(x, ga)
    if supersonic:
        Ma_guess = 1.5 * np.ones_like(A_Acrit)
    else:
        Ma_guess = 0.5 * np.ones_like(A_Acrit)
    return newton(f, Ma_guess, fprime=fp)


def Ma_from_Mash(Mash, ga):
    def f(x):
        return Mash_from_Ma(x, ga) - Mash
    def fp(x):
        return der_Mash_from_Ma(x, ga)
    Ma_guess = 1.5 * np.ones_like(Mash)
    return newton(f, Ma_guess, fprime=fp)


def Ma_from_Posh_Po(Posh_Po, ga):
    def f(x):
        return Posh_Po_from_Ma(x, ga) - Posh_Po
    def fp(x):
        return der_Posh_Po_from_Ma(x, ga)
    Ma_guess = 1.5 * np.ones_like(Posh_Po)
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
    ii = ~np.isinf(Ma)
    V_cpTo = np.sqrt(2.) * np.ones_like(Ma)  # Limit for Ma -> inf
    V_cpTo[ii] = np.sqrt( (ga - 1.0) * Ma[ii] ** 2.
                          / (1. + (ga - 1.) * Ma[ii] **2. /2.) )
    return V_cpTo


def mcpTo_APo_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    mcpTo_APo = (ga / np.sqrt(ga - 1.0) * Ma
                            * To_T ** (-0.5 * (ga + 1.0) / (ga - 1.0)))
    return mcpTo_APo


def mcpTo_AP_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    return ga / np.sqrt(ga - 1.0) * Ma * To_T ** 0.5


def A_Acrit_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    ii = ~(Ma ==.0)  # Do not evaluate when Ma = 0
    A_Acrit = np.ones_like(Ma) * np.inf
    A_Acrit[ii] = ( 1./Ma[ii] * (2. / (ga + 1.0) * To_T[ii])
                    ** (0.5 * (ga + 1.0) / (ga - 1.0)))
    return A_Acrit


def Malimsh_from_ga(ga):
    return np.sqrt((ga - 1.) / ga / 2.) + 0.001


def Mash_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    Malimsh = Malimsh_from_ga(ga)

    # Directly evaluate only when Ma > Malim and not inf
    ii = np.isinf(Ma)
    iv = Ma>Malimsh
    iiv = np.logical_and(~ii, iv)

    Mash = np.ones_like(Ma)  * np.nan  # When Ma < Malim
    Mash[ii] = np.sqrt( (ga - 1.) / ga / 2.) # When Ma -> inf
    Mash[iiv] = (To_T[iiv]
                 / (ga * Ma[iiv] ** 2. - 0.5 * (ga - 1.0))) ** 0.5
    return Mash


def Posh_Po_from_Ma(Ma, ga):
    To_T = To_T_from_Ma(Ma, ga)
    Malimsh = Malimsh_from_ga(ga)

    # Directly evaluate only when Ma > Malim and not inf
    ii = np.isinf(Ma)
    iv = Ma>Malimsh
    iiv = np.logical_and(~ii, iv)

    Posh_Po = np.ones_like(Ma) * np.nan
    Posh_Po[ii] = 0.

    A = 0.5 * (ga + 1.) * Ma[iiv] ** 2. / To_T[iiv]
    B = 2. * ga / (ga + 1.0) * Ma[iiv] ** 2. - 1. / (ga + 1.0) * (ga - 1.0)
    Posh_Po[iiv] = (A ** (ga / (ga - 1.0)) * B ** (-1. / (ga - 1.)))

    return Posh_Po


def from_Ma(var, Ma_in, ga):
    """Evaluate compressible flow quantities as explicit functions of Ma."""

    Ma = np.atleast_1d(Ma_in)
    check_input(Ma, ga)

    # Simple ratios
    if var == 'To_T':
        vout = To_T_from_Ma(Ma, ga)

    elif var == 'Po_P':
        vout = Po_P_from_Ma(Ma, ga)

    elif var == 'rhoo_rho':
        vout = rhoo_rho_from_Ma(Ma, ga)

    # Velocity and mass flow functions
    elif var == 'V_cpTo':
        vout = V_cpTo_from_Ma(Ma, ga)

    elif var == 'mcpTo_APo':
        vout = np.zeros_like(Ma)
        ii = ~np.isinf(Ma)
        vout[ii] = mcpTo_APo_from_Ma(Ma[ii], ga)

    elif var == 'mcpTo_AP':
        vout = mcpTo_AP_from_Ma(Ma, ga)

    # Choking area
    elif var == 'A_Acrit':
        vout = A_Acrit_from_Ma(Ma, ga)

    # Post-shock Mach
    elif var == 'Mash':
        vout = Mash_from_Ma(Ma, ga)

    # Shock pressure ratio
    elif var == 'Posh_Po':
        vout = Posh_Po_from_Ma(Ma, ga)

    else:
        raise ValueError('Incorrect quantity requested')

    if np.size(vout)==1:
        vout = float(vout)

    return vout

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
    dA_Acrit = np.ones_like(Ma) * -np.inf
    ii = ~(Ma == 0.)
    dA_Acrit[ii] = ((2. / (ga + 1.) * To_T[ii]) ** (0.5 * (ga + 1.) / (ga - 1.))
            * (-1./Ma[ii] ** 2. + 0.5 * (ga + 1.) * To_T[ii] ** -1.))
    return dA_Acrit


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
