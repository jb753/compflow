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

* Be as fast as possible, because iterative solutions of
  these equations are often a bottleneck;
* Have a function-oriented interface, with no "ideal gas objects" just to evaluate a formula;
* Operate on non-dimensional quantities only, the one true vocabulary of fluid mechanics.

Benchmarks
**********

