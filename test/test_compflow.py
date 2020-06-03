import compflow as cf
import numpy as np

# Load validation data
ga = 1.4
# Subsonic
var_sub = [
    'Ma',
    'To_T',
    'Po_P',
    'rhoo_rho',
    'V_cpTo',
    'mcpTo_APo',
    'mcpTo_AP',
    'F_mcpTo',
    '4cfL_D',
    'rhoVsq_2Po']

dat_sub = {var_sub[i]: np.loadtxt(
    'test/cued_data_book_subsonic.csv').reshape((-1, 10))[:, i]
    for i in range(len(var_sub))}

# Supersonic
var_sup = var_sub + ['Mash', 'Posh_Po',
                     'Psh_P', 'Posh_P', 'Tsh_T', 'nu', 'Ma2']

dat_sup = {var_sup[i]: np.loadtxt(
    'test/cued_data_book_supersonic.csv').reshape((-1, 17))[3:, i]
    for i in range(len(var_sup))}

# Take reciprocals
for vi in ['To_T', 'Po_P', 'rhoo_rho']:
    dat_sub[vi] = 1. / dat_sub[vi]
    dat_sup[vi] = 1. / dat_sup[vi]

# All variables which the module implements
var_test_sub = ['To_T', 'Po_P', 'rhoo_rho', 'V_cpTo', 'mcpTo_APo', 'mcpTo_AP']
var_test_sup = ['To_T', 'Po_P', 'rhoo_rho', 'V_cpTo', 'mcpTo_APo',
                'mcpTo_AP', 'Mash', 'Posh_Po']

# Number of times to repeat the number crunching
Nrep = 100

#
# Begin test functions

# Check corner cases


def test_Ma_0():
    # Expected values for Ma=0
    val0 = np.array([1., 1., 1., 0., 0., 0., np.nan, np.nan])
    Y0 = {var_test_sup[i]: val0[i] for i in range(len(var_test_sup))}
    for v in var_test_sup:
        # Ma = 0
        assert np.isclose(
            cf.from_Ma(
                v,
                0.,
                ga),
            Y0[v],
            atol=1e-7,
            equal_nan=True)
        if np.isfinite(Y0[v]):
            assert np.isclose(
                cf.to_Ma(
                    v,
                    Y0[v],
                    ga),
                0.0,
                atol=1e-7,
                equal_nan=True)


# Compare code output to CUED Thermofluids Data Book
def test_explicit_sub():
    for v in var_test_sub:
        for _ in range(Nrep):
            Ysub = cf.from_Ma(v, dat_sub['Ma'], ga)
        err_sub = np.abs(Ysub / dat_sub[v] - 1.)
        imax = np.argmax(err_sub)
        assert err_sub[imax] < 0.005, \
            "Ma={0} => {1}={2} with error={3}".format(
                dat_sub['Ma'][imax], v, Ysub[imax], err_sub[imax])


def test_explicit_sup():
    for v in var_test_sup:
        for _ in range(Nrep):
            Ysup = cf.from_Ma(v, dat_sup['Ma'], ga)
        err_sup = np.abs(Ysup / dat_sup[v] - 1.)
        imax = np.argmax(err_sup)
        assert err_sup[imax] < 0.005, \
            "Ma={0} => {1}={2} with error={3}".format(
                dat_sup['Ma'][imax], v, Ysup[imax], err_sup[imax])


def test_inverse_sub():
    for v in var_test_sub:
        for _ in range(Nrep):
            X = cf.to_Ma(v, dat_sub[v], ga)
        err_sub = np.abs(X - dat_sub['Ma'])
        imax = np.argmax(err_sub)
        assert err_sub[imax] <= 0.01, \
            "{0}={1} =>  Ma={2} with error={3}".format(
                v, dat_sub[v][imax], X[imax], err_sub[imax])


def test_inverse_sup():
    for v in var_test_sup:
        for _ in range(Nrep):
            X = cf.to_Ma(v, dat_sup[v], ga, supersonic=True)
        err_sup = np.abs(X - dat_sup['Ma'])
        imax = np.argmax(err_sup)
        assert err_sup[imax] <= 0.01, \
            "{0}={1} =>  Ma={2} with error={3}".format(
                v, dat_sup[v][imax], X[imax], err_sup[imax])


# Check explicit derivative against a numerical approximation
def test_derivative():
    for v in var_test_sup:
        Ma = np.linspace(0., 3.)

        # Explicit
        dYdMa = cf.derivative_from_Ma(v, Ma, ga)

        # Numerical approximation
        Y = cf.from_Ma(v, Ma, ga)
        dYdMa_numerical = np.gradient(Y, Ma)

        # Check discrepancy
        err = dYdMa / dYdMa_numerical - 1.
        assert np.max(err[np.isfinite(err)]) < 0.05
