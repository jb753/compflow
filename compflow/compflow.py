"""Perfect gas compressible flow relations.

This module contains functions to convert back and forth between Mach number
and various other non-dimensional flow quantities.
"""

import numpy as np

from . import native
from .fortran import *

# Initialise empty module-level cache for lookup tables
cache = {}


def _generate_lookup(var, ga, atol):
    """Generate a lookup table for faster inversions to Mach number.

    Args:
        var (str): Name of flow quantity to be used as independent variable.
        ga (float): Ratio of specific heats.
        atol (float): Absolute tolerance on Mach number, default 1e-7.

    Returns:
        f (InterpolatedUnivariateSpline): Callable interpolator that returns
            Mach number as a function of the requested flow quantity."""

    from scipy.interpolate import UnivariateSpline

    # Pick lower limit of table to avoid undefined values
    if var == 'A_Acrit':
        y0 = atol
    else:
        y0 = 0.

    # Start with a small table, uniformly sampled
    Nk = 10
    y = np.linspace(y0, 1., Nk)
    x = from_Ma(var, y, ga)

    # flip if needed
    if x[-1]<x[0]:
        x = np.flip(x)
        y = np.flip(y)

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
        ierr = np.where(err > atol)[0]
        x = np.insert(x,ierr+1,xm[ierr])
        y = np.insert(y,ierr+1,ym[ierr])

    # Add shape restore to the spline object
    def fs(Ma):
        if np.shape(Ma) == ():
            return np.asscalar(f(Ma))
        else:
            return f(Ma)

    return fs

def _cache_lookup(var, ga, atol):
    """Fetch a lookup table from cache, create if needed."""
    if ga not in cache:
        cache[ga] = {}
    if var not in cache[ga]:
        cache[ga][var] = _generate_lookup(var, ga, atol)
    return cache[ga][var]


def to_Ma(var, var_in, ga, supersonic=False):
    """Invert the Mach number relations, solving iteratively if needed."""

    # Choose variable
    if var == 'To_T':
        Ma_out = Ma_from_To_T(var_in, ga)

    elif var == 'Po_P':
        Ma_out = Ma_from_Po_P(var_in, ga)

    elif var == 'rhoo_rho':
        Ma_out = Ma_from_rhoo_rho(var_in, ga)

    elif var == 'V_cpTo':
        Ma_out = Ma_from_V_cpTo(var_in, ga)

    elif var == 'mcpTo_APo':
        Ma_out = Ma_from_mcpTo_APo(var_in, ga, supersonic)

    elif var == 'mcpTo_AP':
        Ma_out = Ma_from_mcpTo_AP(var_in, ga)

    elif var == 'F_mcpTo':
        Ma_out = Ma_from_F_mcpTo(var_in, ga, supersonic)

    elif var == 'A_Acrit':
        Ma_out = Ma_from_A_Acrit(var_in, ga, supersonic)

    elif var == 'Mash':
        Ma_out = Ma_from_Mash(var_in, ga)

    elif var == 'Posh_Po':
        Ma_out = Ma_from_Posh_Po(var_in, ga)

    else:
        raise ValueError('Bad flow quantity requested.')

    return Ma_out


def from_Ma(var, Ma_in, ga):
    """Evaluate compressible flow quantities as explicit functions of Ma."""

    Ma = np.atleast_1d(Ma_in)

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
        vout = mcpTo_APo_from_Ma(Ma, ga)
        # We handle infinite input data explicitly
        vout[np.isinf(Ma)]=0.0

    elif var == 'mcpTo_AP':
        vout = mcpTo_AP_from_Ma(Ma, ga)

    # Impulse
    elif var == 'F_mcpTo':
        vout = F_mcpTo_from_Ma(Ma, ga)

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


def derivative_from_Ma(var, Ma_in, ga):
    """Evaluate quantity derivatives as explicit functions of Mach number."""

    Ma = np.asarray(Ma_in)

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

    # Impulse
    if var == 'F_mcpTo':
        return der_F_mcpTo_from_Ma(Ma, ga)

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


def lookup_mcpTo_APo(mcpTo_APo, ga, atol=1e-6):
    r"""Invert normalised mass flow relation with lookup table.

    This function is equivalent to :func:`compflow.Ma_from_mcpTo_APo` for
    subsonic inputs, but uses a lookup table rather than iteratively solving.
    The lookup table is :ref:`faster <bench>` for large amounts of input data,
    more than about 100 elements.

    The first call of this function with a given value of gamma will
    automatically prepare the table, distributing sampling points as required
    to yield a uniform small error in Mach number. The table is then cached for
    future calls.

    Parameters
    ----------
    mcpTo_APo : array
        Normalised mass flow, :math:`{\dot{m}\sqrt{c_p T_0}}/{A p_0}`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.
    atol : float
        Tolerance on error in inverted :math:`\Ma` values.

    Returns
    -------
    Ma : array
        Mach number, :math:`\Ma`.

    Raises
    ------
    ValueError
        If Mach number is outside the table range of :math:`0\le\Ma\le 1`.

    """
    return _cache_lookup('mcpTo_APo', ga, atol)(mcpTo_APo)

