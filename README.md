# compflow

This module contains functions to convert back and forth between Mach number
and other non-dimensional flow quantities. It is under development but very
much in a usable state.

## Features
* Evaluation of ten non-dimensional flow quantities as explicit functions of Mach number;
* Iteration with Newton's method to invert the explicit relations and solve for Mach number;
* Creation and caching of lookup tables to speed up inversions;
* Fully-vectorised in both directions;
* Tested against a dead-trees data book.

James Brind
June 2020
