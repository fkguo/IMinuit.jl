# IMinuit

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fkguo.github.io/IMinuit.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fkguo.github.io/IMinuit.jl/dev)
[![Build Status](https://github.com/fkguo/IMinuit.jl/workflows/CI/badge.svg)](https://github.com/fkguo/IMinuit.jl/actions)

Julia wrapper of the Python package [`iminuit`](https://github.com/scikit-hep/iminuit). The `minuit` object in `iminuit` is defined as a mutable struct `Fit`.

Functions defined:

`Minuit(fcn; kwds...)`, `Minuit(fcn, start; kwds...)`

 Wrapper of the `iminuit` function `Minuit`.
 `fcn` is the function to be optimized.
 `start`: an array/tuple of the starting values of the parameters.
 `kwds` is the list of keywrod arguments of `Minuit`. For more information, refer to the `iminuit` manual.

Example:
```
    fcn(x) = x[1]^2 + (x[2]-1)^2
    Minuit(fcn, [1, 0]; name = ["a", "b"], error = 0.1*ones(2), fix_a = true, limit_b = (0, 50) )
```
where the parameters are collected in an array `par` which is the argument of `fcn(par)`. In this case,
one can use external code (e.g., `using ForwardDiff: gradient`) to compute the gradient as `gradfun(par) = gradient(fcn, par)`, and include `grad = gradfun` as a keyword argument.

If `fcn` is defined as `fcn(a, b)`, then the starting values need to be set as `Minuit(fcn, a = 1, b = 0)`.

`migrad, minos, hesse, matrix`: wrappers of `iminuit.Minuit.migrad`, `iminuit.Minuit.minos`, `iminuit.Minuit.hesse`, `iminuit.Minuit.matrix`

`func_argnames(f::Function)`:  Extracting the argument names of a function as an array.

```
    Data(x::T, y::T, err::T) where {T<:Vector{Real}}
    Data(df::DataFrame)
```
Fields: `x, y, err, ndata`. This defines a type for data with three columns:` x, y, err`; `ndata` is the number of data rows.
