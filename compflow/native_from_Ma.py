import numpy as np

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

