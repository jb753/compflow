"""This module wraps the fortran backend."""
import numpy as np

import compflow_fort_from_Ma as fort_from_Ma
import compflow_fort_der_from_Ma as fort_der_from_Ma
import compflow_fort_to_Ma as fort_to_Ma

# If the input is a scalar, then put it back to scalar afterwards
def _restore_scalar(func, args):
	if np.shape(args[0]) == ():
		return(func(*args)[0])
	else:
		return(func(*args))

# Functions from Ma
def To_T_from_Ma(Ma, ga):
    r"""Stagnation temperature ratio as function of Mach number.

    .. math::

        \frac{T_0}{T} = 1 + \frac{\gamma - 1}{2} \Ma^2

    Parameters
    ----------
    Ma : array
        Mach number, :math:`\Ma`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    To_T : array
        Stagnation temperature ratio, :math:`T_0/T`.
    """
    return _restore_scalar(fort_from_Ma.to_t, (Ma,ga))

def Po_P_from_Ma(Ma,ga):
    r"""Stagnation pressure ratio as function of Mach number.

    .. math::

        \frac{p_0}{p} = \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)
		^\frac{\gamma}{\gamma - 1}

    Parameters
    ----------
    Ma : array
        Mach number, :math:`\Ma`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    Po_P : array
        Stagnation pressure ratio, :math:`p_0/p`.
    """
    return _restore_scalar(fort_from_Ma.po_p, (Ma,ga))

def rhoo_rho_from_Ma(Ma,ga):
    r"""Stagnation density ratio as function of Mach number.

    .. math::

        \frac{\rho_0}{\rho} = \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)
		^\frac{1}{\gamma - 1}

    Parameters
    ----------
    Ma : array
        Mach number, :math:`\Ma`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    rhoo_rho : array
        Stagnation density ratio, :math:`\rho_0/\rho`.
    """
    return _restore_scalar(fort_from_Ma.rhoo_rho, (Ma,ga))

def V_cpTo_from_Ma(Ma,ga):
    r"""Normalised velocity as function of Mach number.

    .. math::

        \frac{V}{\sqrt{c_p T_0}} = \sqrt{\gamma -1}\, \Ma \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)
		^{-\frac{1}{2}}

    Parameters
    ----------
    Ma : array
        Mach number, :math:`\Ma`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    V_cpTo : array
        Normalised velocity, :math:`V/\sqrt{c_p T_0}`.
    """
    return _restore_scalar(fort_from_Ma.v_cpto, (Ma,ga))

def mcpTo_APo_from_Ma(Ma,ga):
    r"""Normalised mass flow as function of Mach number.

    .. math::

        \frac{\dot{m}\sqrt{c_p T_0}}{A p_0} = \frac{\gamma}{\sqrt{\gamma -1}}\, \Ma \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)
		^{-\frac{1}{2}\frac{\gamma + 1}{\gamma - 1}}

    Parameters
    ----------
    Ma : array
        Mach number, :math:`\Ma`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    mcpTo_APo : array
        Normalised mass flow, :math:`{\dot{m}\sqrt{c_p T_0}}/{A p_0}`.
    """
    return _restore_scalar(fort_from_Ma.mcpto_apo, (Ma,ga))

def mcpTo_AP_from_Ma(Ma,ga):
    r"""Static normalised mass flow as function of Mach number.

    .. math::

        \frac{\dot{m}\sqrt{c_p T_0}}{A p} = \frac{\gamma}{\sqrt{\gamma -1}}\, \Ma \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)
		^{\frac{1}{2}}

    Parameters
    ----------
    Ma : array
        Mach number, :math:`\Ma`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    mcpTo_AP : array
        Static pressure variant of normalised mass flow, :math:`{\dot{m}\sqrt{c_p T_0}}/{A p}`.
    """
    return _restore_scalar(fort_from_Ma.mcpto_ap, (Ma,ga))

def A_Acrit_from_Ma(Ma,ga):
    r"""Ratio of area to choking area as function of Mach number.

    .. math::

        \frac{A}{A_*}  = \frac{1}{\Ma}\left[\frac{2}{\gamma +1} \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)\right]^{\frac{1}{2}\frac{\gamma + 1}{\gamma - 1}}

    Parameters
    ----------
    Ma : array
        Mach number, :math:`\Ma`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    A_Acrit : array
        Ratio of area to choking area, :math:`A/A_*`.
    """
    return _restore_scalar(fort_from_Ma.a_acrit, (Ma,ga))

def Mash_from_Ma(Ma,ga):
    r"""Post-shock Mach number as function of Mach number.

    .. math::

        \Ma_\mathrm{sh}  = \left(\frac{1 + \frac{\gamma - 1}{2} \Ma^2}{\gamma \Ma^2 - \frac{\gamma - 1}{2}} \right)^{\frac{1}{2}}

    Parameters
    ----------
    Ma : array
        Mach number, :math:`\Ma`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    Mash : array
        Post-shock Mach number, :math:`\Ma_\mathrm{sh}`.
    """
    return _restore_scalar(fort_from_Ma.mash, (Ma,ga))

def Posh_Po_from_Ma(Ma,ga):
    r"""Shock stagnation pressure ratio as function of Mach number.

    .. math::

        p_{0\mathrm{sh}}/p_0  = \left(\frac{\frac{\gamma + 1}{2}\Ma^2}{1 + \frac{\gamma - 1}{2} \Ma^2} \right)^{\frac{\gamma}{\gamma-1}}
        \left( \frac{2 \gamma }{\gamma + 1} \Ma^2 - \frac{\gamma - 1}{\gamma + 1}\right)^{\frac{-1}{\gamma -1}}

    Parameters
    ----------
    Ma : array
        Mach number, :math:`\Ma`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    Posh_Po : array
        Shock stagnation pressure ratio, :math:`p_{0\mathrm{sh}}/p_0`.
    """
    return _restore_scalar(fort_from_Ma.posh_po, (Ma,ga))


# Inversions to Ma

def Ma_from_To_T(To_T, ga):
    r"""Mach number as function of stagnation temperature ratio.

    The inverse of :func:`compflow.To_T_from_Ma`, which permits a direct analytical solution.

    .. math::

        \Ma = \sqrt{\frac{2}{\gamma - 1} \left[\frac{T_0}{T} - 1\right]}

    Returns `NaN` if input data is not physically possible, where :math:`{T_0}/{T}<1`.

    Parameters
    ----------
    To_T : array
        Stagnation temperature ratio, :math:`T_0/T`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    Ma : array
        Mach number, :math:`\Ma`.
    """
    return _restore_scalar(fort_to_Ma.to_t, (To_T, ga))

def Ma_from_Po_P(Po_P, ga):
    r"""Mach number as function of stagnation pressure ratio.

    The inverse of :func:`compflow.Po_P_from_Ma`, which permits a direct analytical solution.

    .. math::

        \Ma = \sqrt{ \frac{2}{\gamma - 1}\left[\left(\frac{p_0}{p}\right)^\frac{\gamma - 1}{\gamma} - 1\right]}

    Returns `NaN` if input data is not physically possible, where :math:`{p_0}/{p}<1`.

    Parameters
    ----------
    Po_P : array
        Stagnation pressure ratio, :math:`p_0/p`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    Ma : array
        Mach number, :math:`\Ma`.
    """
    return _restore_scalar(fort_to_Ma.po_p, (Po_P, ga))

def Ma_from_rhoo_rho(rhoo_rho, ga):
    r"""Mach number as function of stagnation density ratio.

    The inverse of :func:`compflow.rhoo_rho_from_Ma`, which permits a direct analytical solution.

    .. math::

        \Ma = \sqrt{ \frac{2}{\gamma - 1}\left[\left(\frac{\rho_0}{\rho}\right)^{\gamma - 1} - 1\right]}

    Returns `NaN` if input data is not physically possible, where :math:`{\rho_0}/{\rho}<1`.

    Parameters
    ----------
    rhoo_rho : array
        Stagnation density ratio, :math:`\rho_0/\rho`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    Ma : array
        Mach number, :math:`\Ma`.
    """
    return _restore_scalar(fort_to_Ma.rhoo_rho, (rhoo_rho, ga))

def Ma_from_V_cpTo(V_cpTo, ga):
    r"""Mach number as function of normalised velocity.

    Inverse of :func:`compflow.V_cpTo_from_Ma`, which permits a direct analytical solution.

    .. math::

        \Ma = \sqrt{\frac{1}{\gamma-1}\left[\frac{\left(\frac{V}{\sqrt{c_p T_0}}\right)^2}{1 - \frac{1}{2}\left(\frac{V}{\sqrt{c_p T_0}}\right)^2}\right]}

    Returns `NaN` if input data is not physically possible, where :math:`V/\sqrt{c_pT_0} < 0` or :math:`V/\sqrt{c_pT_0} > \sqrt{2}`.

    Parameters
    ----------
    V_cpTo : array
        Normalised velocity, :math:`V/\sqrt{c_p T_0}`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    Ma : array
        Mach number, :math:`\Ma`.
    """
    return _restore_scalar(fort_to_Ma.v_cpto, (V_cpTo, ga))

def Ma_from_mcpTo_APo(mcpTo_APo, ga, sup=False):
    r"""Mach number as function of normalised mass flow.

    The inverse of :func:`compflow.mcpTo_APo_from_Ma`, which at a given value of :math:`{\dot{m}\sqrt{c_p T_0}}/{A p_0}` must be solved
    iteratively for :math:`\Ma` using Newton's method.

    .. math::

        \frac{\dot{m}\sqrt{c_p T_0}}{A p_0} = \frac{\gamma}{\sqrt{\gamma -1}}\, \Ma \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)
		^{-\frac{1}{2}\frac{\gamma + 1}{\gamma - 1}}

    For each :math:`{\dot{m}\sqrt{c_p T_0}}/{A p_0}`, there are two possible
    values of :math:`\Ma`. Return the subsonic solution with :math:`\Ma\le 1`
    by default; the supersonic solution with :math:`\Ma>`` is retrived by
    setting the parameter `sup=True`.

    Returns `NaN` if input data is not physically possible, where
    :math:`{\dot{m}\sqrt{c_p T_0}}/{A p_0} < 0`. The normalised mass flow
    reaches a maximum at the sonic velocity :math:`\Ma=1`. Input data above the
    maximum value correspond to choking --- also return `NaN` in this case.

    Parameters
    ----------
    mcpTo_APo : array
        Normalised mass flow, :math:`{\dot{m}\sqrt{c_p T_0}}/{A p_0}`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.
    sup : bool, default False
        If true, return the supersonic solution, otherwise the subsonic solution.

    Returns
    -------
    Ma : array
        Mach number, :math:`\Ma`.
    """
    return _restore_scalar(fort_to_Ma.mcpto_apo, (mcpTo_APo, ga, sup))

def Ma_from_mcpTo_AP(mcpTo_AP, ga):
    r"""Mach number as function of static normalised mass flow.

    The inverse of :func:`compflow.mcpTo_AP_from_Ma`, which at a given value of
    :math:`{\dot{m}\sqrt{c_p T_0}}/{A p}` must be solved
    iteratively for :math:`\Ma` using Newton's method.

    .. math::

        \frac{\dot{m}\sqrt{c_p T_0}}{A p} = \frac{\gamma}{\sqrt{\gamma -1}}\, \Ma \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)
		^{\frac{1}{2}}

    Returns `NaN` if input data is not physically possible, where
    :math:`{\dot{m}\sqrt{c_p T_0}}/{A p} < 0`.

    Parameters
    ----------
    mcpTo_AP : array
        Static normalised mass flow, :math:`{\dot{m}\sqrt{c_p T_0}}/{A p}`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    Ma : array
        Mach number, :math:`\Ma`.
    """
    return _restore_scalar(fort_to_Ma.mcpto_ap, (mcpTo_AP, ga))

def Ma_from_A_Acrit(A_Acrit, ga):
    r"""Mach number as function of area to choking area ratio.

    The inverse of :func:`compflow.A_Acrit_from_Ma`, which at a given value of
    :math:`A/A_*` must be solved iteratively for :math:`\Ma` using Newton's
    method.

    .. math::

        \frac{A}{A_*}  = \frac{1}{\Ma}\left[\frac{2}{\gamma +1} \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)\right]^{\frac{1}{2}\frac{\gamma + 1}{\gamma - 1}}

    Returns `NaN` if input data is not physically possible, where
    :math:`A/A_* < 1`.

    Parameters
    ----------
    A_Acrit : array
        Ratio of area to choking area, :math:`A/A_*`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    Ma : array
        Mach number, :math:`\Ma`.
    """
    return _restore_scalar(fort_to_Ma.a_acrit, (A_Acrit, ga))

def Ma_from_Mash(Mash, ga):
    r"""Mach number as function of post-shock Mach number.

    The inverse of :func:`compflow.Mash_from_Ma`, which at a given value of
    :math:`\Ma_\mathrm{sh}` must be solved iteratively for :math:`\Ma` using
    Newton's method.

    .. math::

        \Ma_\mathrm{sh}  = \left(\frac{1 + \frac{\gamma - 1}{2} \Ma^2}{\gamma \Ma^2 - \frac{\gamma - 1}{2}} \right)^{\frac{1}{2}}

    Returns `NaN` if input data is not physically possible, where
    :math:`\Ma_\mathrm{sh}>1` or :math:`\Ma_\mathrm{sh}<0`.

    Parameters
    ----------
    Mash : array
        Post-shock Mach number, :math:`\Ma_\mathrm{sh}`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    Ma : array
        Mach number, :math:`\Ma`.
    """
    return _restore_scalar(fort_to_Ma.mash, (Mash, ga))

def Ma_from_Posh_Po(Posh_Po, ga):
    r"""Mach number as function of shock stagnation pressure ratio.

    The inverse of :func:`compflow.Posh_Po_from_Ma`, which at a given value of
    :math:`p_{0\mathrm{sh}}/p_0` must be solved iteratively for :math:`\Ma`
    using Newton's method.

    .. math::

        \frac{p_{0\mathrm{sh}}}{p_0}  = \left(\frac{\frac{\gamma + 1}{2}\Ma^2}{1 + \frac{\gamma - 1}{2} \Ma^2} \right)^{\frac{\gamma}{\gamma-1}}
        \left( \frac{2 \gamma }{\gamma + 1} \Ma^2 - \frac{\gamma - 1}{\gamma + 1}\right)^{\frac{-1}{\gamma -1}}

    Returns `NaN` if input data is not physically possible, where
    :math:`p_{0\mathrm{sh}}/p_0 > 1` or :math:`p_{0\mathrm{sh}}/p_0 < 0`.

    Parameters
    ----------
    Posh_Po : array
        Shock stagnation pressure ratio, :math:`p_{0\mathrm{sh}}/p_0`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    Ma : array
        Mach number, :math:`\Ma`.
    """
    return _restore_scalar(fort_to_Ma.posh_po, (Posh_Po, ga))

# Derivatives from Ma
der_To_T_from_Ma = fort_der_from_Ma.to_t
der_Po_P_from_Ma = fort_der_from_Ma.po_p
der_rhoo_rho_from_Ma = fort_der_from_Ma.rhoo_rho
der_V_cpTo_from_Ma = fort_der_from_Ma.v_cpto
der_mcpTo_APo_from_Ma = fort_der_from_Ma.mcpto_apo
der_mcpTo_AP_from_Ma = fort_der_from_Ma.mcpto_ap
der_A_Acrit_from_Ma = fort_der_from_Ma.a_acrit
der_Mash_from_Ma = fort_der_from_Ma.mash
der_Posh_Po_from_Ma = fort_der_from_Ma.posh_po
