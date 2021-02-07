"""This module wraps the fortran backend."""

import compflow_fort_from_Ma as fort_from_Ma
import compflow_fort_der_from_Ma as fort_der_from_Ma
import compflow_fort_to_Ma as fort_to_Ma

# Functions from Ma
def To_T_from_Ma(Ma, ga):
    r"""Stagnation temperature ratio as function of Mach Number.

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
    return fort_from_Ma.to_t(Ma,ga)

def Po_P_from_Ma(Ma,ga):
    r"""Stagnation pressure ratio as function of Mach Number.

    .. math::

        \frac{p_0}{p} = \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)
		^\dfrac{\gamma}{\gamma - 1}

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
    return fort_from_Ma.po_p(Ma,ga)

def rhoo_rho_from_Ma(Ma,ga):
    r"""Stagnation density ratio as function of Mach Number.

    .. math::

        \frac{\rho_0}{\rho} = \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)
		^\dfrac{1}{\gamma - 1}

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
    return fort_from_Ma.rhoo_rho(Ma,ga)

def V_cpTo_from_Ma(Ma,ga):
    r"""Normalised velocity as function of Mach Number.

    .. math::

        \frac{V}{\sqrt{c_p T_0}} = \sqrt{\gamma -1}\, \Ma \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)
		^{-\dfrac{1}{2}}

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
    return fort_from_Ma.v_cpto(Ma,ga)


def mcpTo_APo_from_Ma(Ma,ga):
    r"""Normalised mass flow as function of Mach number.

    .. math::

        \frac{\dot{m}\sqrt{c_p T_0}}{A p_0} = \frac{\gamma}{\sqrt{\gamma -1}}\, \Ma \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)
		^{-\dfrac{1}{2}\dfrac{\gamma + 1}{\gamma - 1}}

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
    return fort_from_Ma.mcpto_apo(Ma,ga)

def mcpTo_AP_from_Ma(Ma,ga):
    r"""Static pressure variant of normalised mass flow as function of Mach number.

    .. math::

        \frac{\dot{m}\sqrt{c_p T_0}}{A p} = \frac{\gamma}{\sqrt{\gamma -1}}\, \Ma \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)
		^{\dfrac{1}{2}}

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
    return fort_from_Ma.mcpto_ap(Ma,ga)

def A_Acrit_from_Ma(Ma,ga):
    r"""Ratio of area to choking area as function of Mach number.

    .. math::

        \frac{A}{A_*}  = \frac{1}{\Ma}\left[\frac{2}{\gamma +1} \left(1 + \frac{\gamma - 1}{2} \Ma^2 \right)\right]^{\dfrac{1}{2}\dfrac{\gamma + 1}{\gamma - 1}}

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
    return fort_from_Ma.a_acrit(Ma,ga)

def Mash_from_Ma(Ma,ga):
    r"""Post-shock Mach number as function of (pre-shock) Mach number.

    .. math::

        \Ma_\mathrm{sh}  = \left(\frac{1 + \dfrac{\gamma - 1}{2} \Ma^2}{\gamma \Ma^2 - \dfrac{\gamma - 1}{2}} \right)^{\dfrac{1}{2}}

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
    return fort_from_Ma.mash(Ma,ga)


def Posh_Po_from_Ma(Ma,ga):
    r"""Shock stagnation pressure ratio as function of (pre-shock) Mach number.

    .. math::

        p_{0\mathrm{sh}}  = \left(\frac{\dfrac{\gamma + 1}{2}\Ma^2}{1 + \dfrac{\gamma - 1}{2} \Ma^2} \right)^{\dfrac{\gamma}{\gamma-1}}
        \left( \frac{2 \gamma }{\gamma + 1} \Ma^2 - \frac{\gamma - 1}{\gamma + 1}\right)^{\dfrac{-1}{\gamma -1}}

    Parameters
    ----------
    Ma : array
        Mach number, :math:`\Ma`.
    ga : float
        Ratio of specific heats, :math:`\gamma`.

    Returns
    -------
    Mash : array
        Shock stagnation pressure ratio, :math:`p_{0\mathrm{sh}}/p_0`.
    """
    return fort_from_Ma.posh_po(Ma,ga)

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

# Inversions to Ma
Ma_from_To_T = fort_to_Ma.to_t
Ma_from_Po_P = fort_to_Ma.po_p
Ma_from_rhoo_rho = fort_to_Ma.rhoo_rho
Ma_from_V_cpTo = fort_to_Ma.v_cpto
Ma_from_mcpTo_APo = fort_to_Ma.mcpto_apo
Ma_from_mcpTo_AP = fort_to_Ma.mcpto_ap
Ma_from_A_Acrit = fort_to_Ma.a_acrit
Ma_from_Mash = fort_to_Ma.mash
Ma_from_Posh_Po = fort_to_Ma.posh_po
