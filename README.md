# compflow

The compflow library contains functions to convert back and forth between Mach
number and other non-dimensional groups in compressible flows. By using a
NumPy--Fortran interface, the code is vectorised and lightning-fast, yielding a
speed-up of up to two orders of magnitude.

[Full documentation](https://jamesbrind.uk/compflow-docs/) is available online. 

![Compressible flow quantities](sample/sample.png)

## Features

* Evaluation of ten non-dimensional flow quantities as explicit functions of
  Mach number;
* Iteration with Newton's method to invert explicit relations and solve for
  Mach number;
* Creation and caching of lookup tables to speed up inversions;
* Fortran-accelerated, fully-vectorised in both directions.

## Basic usage
 
compflow is available on the Python Package Index, so installation is as simple as,

```python

   python3 -m pip install compflow

   ```
**Note:** as the library uses Fortran subroutines behind the scenes, you will
need a working Fortran compiler for the installation to complete successfully. 

We can now start doing some calculations. First, an explicit evaluation of
stagnation pressure ratio given a Mach number,

```python

   >>> import compflow
   >>> ga = 1.4
   >>> compflow.Po_P_from_Ma(0.3, ga)
   1.0644302861529382

```

Second, an inversion of flow function where
iterative solution for Mach number is required,

```python

   >>> compflow.Ma_from_mcpTo_APo(0.8, ga)
   0.39659360325173604

```

The names and symbols of non-dimensional quantities are fairly
self-explanatory, but a full list is given in the
[Nomenclature](https://jamesbrind.uk/compflow-docs/api.html#nomenclature).
All functions and the equations used for the calculations are documented in the
[API](https://jamesbrind.uk/compflow-docs/api.html).

Numpy arrays are also accepted as inputs,

```python

   >>> import numpy
   >>> Ma1 = numpy.array([0., 0.5, 1., 2.])
   >>> compflow.To_T_from_Ma(Ma1, ga)
   array([1.  , 1.05, 1.2 , 1.8 ])
   >>> Ma2 = numpy.array([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]])
   >>> compflow.To_T_from_Ma(Ma2, ga)
   array([[1.002, 1.008],
          [1.018, 1.032],
          [1.05 , 1.072]])
```

When solving for Mach number at a given normalised mass flow, it is assumed
that we are on the subsonic branch of the curve unless a flag is specified.
Where no solution is possible, i.e. if the flow would choke, `NaN` is
returned,

```python

   >>> capacity = [0.6, 2.]
   >>> compflow.Ma_from_mcpTo_APo(capacity, ga)
   array([0.28442265,        nan])
   >>> compflow.Ma_from_mcpTo_APo(capacity, ga, sup=True)
   array([2.27028708,        nan])
```

## TODO

* Sort out packaging so that NumPy gets installed automatically (distutils due to be deprecated?).

James Brind
Mar 2022
