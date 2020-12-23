import numpy as np
from scipy.optimize import newton
from .native_from_Ma import *

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

