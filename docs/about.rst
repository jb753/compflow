About
=====

Motivation
**********

compflow was written because the existing options did not meet the author's
needs. In particular, for internal flow problems, a very common operation is
solving for Mach number given a normalised mass flow
:math:`\dot{m}\sqrt{c_pT_0}/Ap_0` (also known as flow function or capacity).
Surprisingly, although there are many similar aerodynamics libraries available,
none had this feature.

Other distinctive aspects of the design philosophy are that the library should:

* Be as fast as possible, because iterative solutions of these equations are
  often a bottleneck;
* Operate on vectors of input data, for scalability;
* Have a function-oriented interface, with no "ideal gas objects" just to
  evaluate a formula;
* Operate on non-dimensional quantities only, the one true vocabulary of fluid
  mechanics.

.. _bench:

Benchmarks
**********

This section presents benchmarks that quantify the speed-up of compflow
functions, written in Fortran, compared to a native NumPy implementation. The
input data can be vectorised, so the benchmarks is repeated for different sized
arrays. The script which generates these benchmarks is `available on the GitHub
<https://github.com/jb753/compflow/blob/master/test/run_bench.py>`_.

.. image:: _static/bench_forward.svg
   :align: center
   :width: 60%

Considering forward evaluation of normalised mass flow,
:func:`compflow.mcpTo_APo_from_Ma` delivers a speed-up of 5.4 for scalar
inputs, which is worth having. The native implementation requires multiple
calls to NumPy functions. Although each NumPy function is individually
well-optimised C, there is an advantage to doing the entire calculation in one
compiled go. For large arrays, the relative overhead of going in and out of
compiled functions decreases, and the speed-up diminishes until the
implementations are equal.

The normalised mass flow relation does not admit an analytical solution for
Mach number, so either iteration is required, or a lookup table can be stored
for many values of :math:`\dot{m}\sqrt{c_pT_0}/Ap_0`. The library function
:func:`compflow.Ma_from_mcpTo_APo` uses iteration, but compflow also implements
lookup tables using SciPy splines. The native implementation uses the SciPy
Newton solver.

.. image:: _static/bench_inverse.svg
   :align: center
   :width: 60%

For scalar inputs, the lookup table yields a speed-up of 46 compared to the
native implementation. This does not include the computations required to
initialise the lookup table. The compiled Fortran version from compflow yields
a larger speed-up of 180, despite the fact that it re-calculates the values
every time! The advantage comes from the specialised nature of the compflow
code --- the SciPy spline and Newton solver are complex, general tools with
quite a bit of extra logic that is not used here. However, for large arrays of
input data, the lookup table scales better, with a speed-up of 11 compared to
3.2 for the Fortran implementation.
