Basic Usage
===========

compflow is available on the Python Package Index, so installation is as simple as,

.. code-block:: bash

   python3 -m pip install compflow

.. note::

   As the library uses Fortran subroutines behind the scenes, you will need a
   working Fortran compiler for the installation to complete successfully. 

.. testsetup:: *

   import compflow

We can now start doing some calculations. First, an explicit evaluation of
stagnation pressure ratio :math:`p_0/p` given a Mach number :math:`\Ma`,

.. doctest::

   >>> import compflow
   >>> ga = 1.4
   >>> compflow.Po_P_from_Ma(0.3, ga)
   1.0644302861529382

Second, an inversion of flow function :math:`\dot{m}\sqrt{c_pT_0}/Ap_0` where
iterative solution for :math:`\Ma` is required,

.. doctest::

   >>> compflow.Ma_from_mcpTo_APo(0.8, ga)
   0.39659360325173604

The names and symbols of non-dimensional quantities are fairly
self-explanatory, but a full list is given in the :ref:`nomen`. All functions
and the equations used for the calculations are documented in the :ref:`api`.

Numpy arrays are also accepted as inputs,

.. doctest::

   >>> import numpy
   >>> Ma = numpy.array([0., 0.5, 1., 2.])
   >>> compflow.To_T_from_Ma(Ma, ga)
   array([1.  , 1.05, 1.2 , 1.8 ])

When solving for Mach number at a given normalised mass flow, it is assumed
that we are on the subsonic branch of the curve unless a flag is specified.
Where no solution is possible, for example if the flow would choke, `NaN` is
returned,

.. doctest::

   >>> capacity = [0.6, 2.]
   >>> compflow.Ma_from_mcpTo_APo(capacity, ga)
   array([0.28442265,        nan])
   >>> compflow.Ma_from_mcpTo_APo(capacity, ga, sup=True)
   array([2.27028708,        nan])
