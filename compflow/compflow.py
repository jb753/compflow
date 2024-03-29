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
    if var == "A_Acrit":
        y0 = atol
    else:
        y0 = 0.0

    # Start with a small table, uniformly sampled
    Nk = 10
    y = np.linspace(y0, 1.0, Nk)
    x = from_Ma(var, y, ga)

    # flip if needed
    if x[-1] < x[0]:
        x = np.flip(x)
        y = np.flip(y)

    # Loop until we reach tolerance
    err_max = np.inf
    N_max = 50
    for n in range(N_max):

        # Make interpolator
        f = UnivariateSpline(x, y, k=1, s=0.0, ext="raise")

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
        x = np.insert(x, ierr + 1, xm[ierr])
        y = np.insert(y, ierr + 1, ym[ierr])

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
    if var == "To_T":
        Ma_out = Ma_from_To_T(var_in, ga)

    elif var == "Po_P":
        Ma_out = Ma_from_Po_P(var_in, ga)

    elif var == "rhoo_rho":
        Ma_out = Ma_from_rhoo_rho(var_in, ga)

    elif var == "V_cpTo":
        Ma_out = Ma_from_V_cpTo(var_in, ga)

    elif var == "mcpTo_APo":
        Ma_out = Ma_from_mcpTo_APo(var_in, ga, supersonic)

    elif var == "mcpTo_AP":
        Ma_out = Ma_from_mcpTo_AP(var_in, ga)

    elif var == "F_mcpTo":
        Ma_out = Ma_from_F_mcpTo(var_in, ga, supersonic)

    elif var == "A_Acrit":
        Ma_out = Ma_from_A_Acrit(var_in, ga, supersonic)

    elif var == "Mash":
        Ma_out = Ma_from_Mash(var_in, ga)

    elif var == "Posh_Po":
        Ma_out = Ma_from_Posh_Po(var_in, ga)

    else:
        raise ValueError("Bad flow quantity requested.")

    return Ma_out


def from_Ma(var, Ma_in, ga):
    """Evaluate compressible flow quantities as explicit functions of Ma."""

    Ma = np.atleast_1d(Ma_in)

    # Simple ratios
    if var == "To_T":
        vout = To_T_from_Ma(Ma, ga)

    elif var == "Po_P":
        vout = Po_P_from_Ma(Ma, ga)

    elif var == "rhoo_rho":
        vout = rhoo_rho_from_Ma(Ma, ga)

    # Velocity and mass flow functions
    elif var == "V_cpTo":
        vout = V_cpTo_from_Ma(Ma, ga)

    elif var == "mcpTo_APo":
        vout = mcpTo_APo_from_Ma(Ma, ga)
        # We handle infinite input data explicitly
        vout[np.isinf(Ma)] = 0.0

    elif var == "mcpTo_AP":
        vout = mcpTo_AP_from_Ma(Ma, ga)

    # Impulse
    elif var == "F_mcpTo":
        vout = F_mcpTo_from_Ma(Ma, ga)

    # Choking area
    elif var == "A_Acrit":
        vout = A_Acrit_from_Ma(Ma, ga)

    # Post-shock Mach
    elif var == "Mash":
        vout = Mash_from_Ma(Ma, ga)

    # Shock pressure ratio
    elif var == "Posh_Po":
        vout = Posh_Po_from_Ma(Ma, ga)

    else:
        raise ValueError("Incorrect quantity requested")

    if np.size(vout) == 1:
        vout = float(vout)

    return vout


def derivative_from_Ma(var, Ma_in, ga):
    """Evaluate quantity derivatives as explicit functions of Mach number."""

    Ma = np.asarray(Ma_in)

    # Simple ratios
    if var == "To_T":
        return der_To_T_from_Ma(Ma, ga)

    if var == "Po_P":
        return der_Po_P_from_Ma(Ma, ga)

    if var == "rhoo_rho":
        return der_rhoo_rho_from_Ma(Ma, ga)

    # Velocity and mass flow functions
    if var == "V_cpTo":
        return der_V_cpTo_from_Ma(Ma, ga)

    if var == "mcpTo_APo":
        return der_mcpTo_APo_from_Ma(Ma, ga)

    if var == "mcpTo_AP":
        return der_mcpTo_AP_from_Ma(Ma, ga)

    # Impulse
    if var == "F_mcpTo":
        return der_F_mcpTo_from_Ma(Ma, ga)

    # Choking area
    if var == "A_Acrit":
        return der_A_Acrit_from_Ma(Ma, ga)

    # Post-shock Mack number
    if var == "Mash":
        return der_Mash_from_Ma(Ma, ga)

    # Shock pressure ratio
    if var == "Posh_Po":
        return der_Posh_Po_from_Ma(Ma, ga)

    # Throw an error if we don't recognise the requested variable
    raise ValueError("Invalid quantity requested: {}.".format(var))


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
    return _cache_lookup("mcpTo_APo", ga, atol)(mcpTo_APo)


def static_from_stagnation(To, Po, Ma, ga):
    r"""Given dimensional stagnation state and Mach number, get static state.

    A convenience function that makes calls to :func:`compflow.To_T_from_Ma`
    and :func:`compflow.Po_P_from_Ma` in order to convert a dimensional
    stagnation pressure and temperature to static quantities (in the same
    units).

    This function is the inverse of :func:`compflow.stagnation_from_static` and
    a special case of :func:`compflow.change_frame`.

    Parameters
    ----------
    To : array
        Stagnation temperature, :math:`T_0`.
    Po : array
        Stagnation pressure, :math:`p_0`.
    Ma : array
        Mach number, :math:`\Ma`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    T : array
        Static temperature, :math:`T`.
    P : array
        Static pressure, :math:`p`.

    """

    return change_frame(To, Po, Ma, 0.0, ga)


def stagnation_from_static(T, P, Ma, ga):
    r"""Given dimensional static state and Mach number, get stagnation state.

    A convenience function that makes calls to :func:`compflow.To_T_from_Ma`
    and :func:`compflow.Po_P_from_Ma` in order to convert a dimensional
    stagnation pressure and temperature to static quantities (in the same
    units).

    This function is the inverse of :func:`compflow.static_from_stagnation` and
    a special case of :func:`compflow.change_frame`.

    Parameters
    ----------
    T : array
        Static temperature, :math:`T`.
    P : array
        Static pressure, :math:`p`.
    Ma : array
        Mach number, :math:`\Ma`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    To : array
        Stagnation temperature, :math:`T_0`.
    Po : array
        Stagnation pressure, :math:`p_0`.

    """

    return change_frame(T, P, 0.0, Ma, ga)


def change_frame(To1, Po1, Ma1, Ma2, ga):
    r"""Convert a dimensional stagnation state between reference frames.

    Suppose we have a stagnation state :math:`T_{01}, P_{01}` in one reference
    frame, with a known Mach number relative to that frame :math:`\Ma_1`. Then,
    given the Mach number relative to a *different* reference frame
    :math:`\Ma_2`, we would like to know the stagnation state in the new
    reference frame :math:`T_{02}, P_{02}`.

    The static state is independent of reference frame. So we convert the
    stagnation state to static and back again in the new reference frame,
    calling :func:`compflow.To_T_from_Ma` and :func:`compflow.Po_P_from_Ma`
    with the respective Mach numbers.

    Swapping the input Mach numbers will invert this function.

    Parameters
    ----------
    To1 : array
        Stagnation temperature in frame 1, :math:`T_{01}`.
    Po1 : array
        Stagnation pressure in frame 1, :math:`p_{01}`.
    Ma1 : array
        Mach number relative to frame 1, :math:`\Ma_1`.
    Ma2 : array
        Mach number relative to frame 2, :math:`\Ma_2`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    To2 : array
        Stagnation temperature in frame 2, :math:`T_{02}`.
    Po2 : array
        Stagnation pressure in frame 2, :math:`p_{02}`.
    """

    # Get ratios
    Po1_P = Po_P_from_Ma(Ma1, ga)
    To1_T = To_T_from_Ma(Ma1, ga)
    Po2_P = Po_P_from_Ma(Ma2, ga)
    To2_T = To_T_from_Ma(Ma2, ga)

    # Convert
    To2 = To1 / To1_T * To2_T
    Po2 = Po1 / Po1_P * Po2_P

    return To2, Po2


def specific_heats(ga, rgas):
    r"""Specific heats from gas constant and specific heat ratio.

    For a perfect gas, the specific heat capacities at constant pressure and
    volume are given by,

    .. math::

         c_p  = R \frac{\gamma}{\gamma - 1} \, , \quad  c_v  = R \frac{1}{\gamma - 1}

    Parameters
    ----------
    ga : float
        Ratio of specific heats, :math:`\gamma`.
    rgas : float
        Specific gas constant, :math:`R`.

    Returns
    -------
    cp : float
        Specific heat capacity at constant pressure, :math:`c_p`.
    cv : float
        Specific heat capacity at constant volume, :math:`c_v`.
    """
    cv = rgas / (ga - 1.0)
    cp = cv * ga
    return cp, cv


def gas_constant(ga, cp):
    r"""Specific gas constant from specific heat ratio and heat capacity.

    For a perfect gas, the specific gas constant is,

    .. math::

        R = c_p \frac{\gamma - 1}{\gamma}

    Parameters
    ----------
    ga : float
        Ratio of specific heats, :math:`\gamma`.
    cp : float
        Specific heat capacity at constant pressure, :math:`c_p`.

    Returns
    -------
    rgas : float
        Specific gas constant, :math:`R`.
    """
    return (ga - 1.0) / ga * cp


def primitive_from_conserved(r, ro, rovx, rovr, rorvt, roe, ga, rgas):
    r"""Get velocity, pressure temperature from conserved quantities.

    Computational fluid dynamics codes often store flow fields in terms of
    *conserved* quantities, the net fluxes of which sum to zero for each cell:

    .. math::

        \rho \, ,\ \rho V_x \, ,\ \rho V_r \, ,\ \rho r V_\theta \, ,\ \rho e \, ,

    where the specific internal energy :math:`e = c_v T + \tfrac{1}{2} V^2`.
    Of interest for further analysis are an equivalent set of five *primitive*
    quantities that are easier to work with:

    .. math::

        V_x \, ,\ V_r \, ,\ V_\theta \, ,\ p \, ,\ T \, .

    This function converts the *conserved* quantities to *primitive* quantities
    by dividing out density and using the ideal gas law. The inverse of this
    function is :func:`compflow.conserved_from_primitive`.

    Parameters
    ----------
    r : array
        Radial coordinate, :math:`r`.
    ro : array
        Density, :math:`\rho`.
    rovx : array
        Axial momentum flux, :math:`\rho V_x`.
    rovr : array
        Radial momentum flux, :math:`\rho V_r`.
    rorvt : array
        Moment of angular momentum flux, :math:`\rho r V_\theta`.
    rovx : array
        Volumetric internal energy, :math:`\rho e`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.
    rgas : float
        Specific gas constant, :math:`R`.

    Returns
    -------
    vx : array
        Axial velocity, :math:`V_x`.
    vr : array
        Radial velocity, :math:`V_r`.
    vt : array
        Circumferential velocity, :math:`V_\theta`.
    P : array
        Static pressure, :math:`p`.
    T : array
        Static temperature, :math:`T`.
    """

    cp, cv = specific_heats(ga, rgas)

    # Divide out density
    vx = rovx / ro
    vr = rovr / ro
    vt = rorvt / ro / r
    e = roe / ro

    # Calculate secondary variables
    vsq = vx ** 2.0 + vr ** 2.0 + vt ** 2.0
    T = (e - 0.5 * vsq) / cv
    P = ro * rgas * T

    return vx, vr, vt, P, T


def conserved_from_primitive(r, vx, vr, vt, P, T, ga, rgas):
    r"""Get conserved quantities from velocity, pressure, temperature.

    Two thermodynamic properties and three velocity components are sufficient
    to specify a three-dimensional flow field. We will call these the
    *primitive* set of quantities,

    .. math::

        V_x \, ,\ V_r \, ,\ V_\theta \, ,\ p \, ,\ T \, .

    However, computational fluid dynamics codes often store flow fields in
    terms of *conserved* quantities, the net fluxes of which sum to zero for
    each cell:

    .. math::

        \rho \, ,\ \rho V_x \, ,\ \rho V_r \, ,\ \rho r V_\theta \, ,\ \rho e \, ,

    where the specific internal energy :math:`e = c_v T + \tfrac{1}{2} V^2`.

    This function converts *primitive* quantities to *conserved* quantities,
    and is the inverse of :func:`compflow.primitive_from_conserved`.

    Parameters
    ----------
    r : array
        Radial coordinate, :math:`r`.
    vx : array
        Axial velocity, :math:`V_x`.
    vr : array
        Radial velocity, :math:`V_r`.
    vt : array
        Circumferential velocity, :math:`V_\theta`.
    P : array
        Static pressure, :math:`p`.
    T : array
        Static temperature, :math:`T`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.
    rgas : float
        Specific gas constant, :math:`R`.

    Returns
    -------
    ro : array
        Density, :math:`\rho`.
    rovx : array
        Axial momentum flux, :math:`\rho V_x`.
    rovr : array
        Radial momentum flux, :math:`\rho V_r`.
    rorvt : array
        Moment of angular momentum flux, :math:`\rho r V_\theta`.
    roe : array
        Volumetric internal energy, :math:`\rho e`.
    """

    cp, cv = specific_heats(ga, rgas)
    vsq = vx ** 2.0 + vr ** 2.0 + vt ** 2.0
    ro = P / rgas / T
    rovx = ro * vx
    rovr = ro * vr
    rorvt = ro * r * vt
    roe = ro * (cv * T + 0.5 * vsq)
    return ro, rovx, rovr, rorvt, roe
