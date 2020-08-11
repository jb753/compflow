# compflow

This library contains functions to convert back and forth between Mach number
and other non-dimensional compressible flow quantities.

![Compressible flow quantities](sample/sample.svg)

The features are more-or-less complete and development is now focussed on
improving performance. The API should remain stable, but an absence of breakage
in future versions is not guaranteed.

## Features

* Evaluation of ten non-dimensional flow quantities as explicit functions of
  Mach number;
* Iteration with Newton's method to invert the explicit relations and solve for
  Mach number;
* Creation and caching of lookup tables to speed up inversions;
* Fully-vectorised in both directions, accepts NumPy arrays.

## Usage

compflow is available on the Python Package Index, so installation is as simple as,
```
   $ python3 -m pip install compflow
```

We can then start doing some calculations,
```python
>>> import compflow
>>> compflow.from_Ma(var='Po_P', Ma_in=0.3, ga=1.4)  # Explicit evaluation
1.0644302861529382
>>> compflow.to_Ma(var='Po_P', var_in=1.064, ga=1.4)  # Iterative solve
0.2990184500710573
```
The functions `from_Ma` and `to_Ma` are the primary interface to the library.
They take a string argument indicating which quantity should be found from or
converted to Mach number; see the docs for valid choices.

Numpy arrays are also accepted as inputs,
```python
>>> import numpy
>>> Ma = numpy.array([0., 0.5, 1., 2.])
>>> compflow.from_Ma(var='Po_P', Ma_in=Ma, ga=1.4)
array([1.        , 1.18621264, 1.89292916, 7.82444907])
```

When solving for Mach number, it is assumed that we are on the subsonic branch
of the curve unless a flag is specified. Where no solution is possible, for
example if the flow would choke, `nan` is returned,
```python 
>>> capacity = [0.6, 2.]
>>> compflow.to_Ma(var='mcpTo_APo', var_in=capacity, ga=1.4)
array([0.28442265,        nan])
>>> compflow.to_Ma(var='mcpTo_APo', var_in=capacity, ga=1.4, supersonic=True)
array([2.27028708,        nan])
```

Lower level functions for each quantity are also available, but these are fussy
about their inputs: some accept only Numpy arrays.
```python
>>> compflow.Po_P_from_Ma(0.3,1.4)
1.064430286152938
>>> compflow.A_Acrit_from_Ma(numpy.atleast_1d(0.3),1.4)
array([2.03506526])
```

## Why compflow?

Like many libraries, compflow was written because the existing options did not
meet the authors' needs.

* Implements non-dimensional mass flow or capacity, which is commonly used for
  internal flow calculations;
* Simple, function-oriented interface;
* Supports inversions to solve for Mach number;
* Fully non-dimensional, no units-related nonsense.

## TODO

* Investigate speedups using `numexpr`, `Numba`, and `PyPy`

James Brind
August 2020
