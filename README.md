# IMinuit

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fkguo.github.io/IMinuit.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fkguo.github.io/IMinuit.jl/dev)
[![Build Status](https://github.com/fkguo/IMinuit.jl/workflows/CI/badge.svg)](https://github.com/fkguo/IMinuit.jl/actions)

Julia wrapper of the Python package [`iminuit`](https://github.com/scikit-hep/iminuit), which is the interface to the C++ MINUIT2, widely used in fitting in the high-energy physics community. The `minuit` object in `iminuit` is defined as an `AbstractFit`:
if using array parameters, then `ArrayFit`;
if using individual parameters, then `Fit`.

Install by `]add https://github.com/fkguo/IMinuit.jl`

Functions defined:

## Functions in `iminuit`


```
Minuit(fcn; kwds...)::Fit
Minuit(fcn, start; kwds...)::ArrayFit
```

 Wrapper of the `iminuit` function `Minuit`.
 `fcn` is the function to be optimized.
 `start`: an array/tuple of the starting values of the parameters.
 `kwds` is the list of keywrod arguments of `Minuit`. For more information, refer to the `iminuit` manual.
 The `Fit` one, for which `fcn` takes individual parameters as variables,  is generally
 faster than the `ArrayFit` one, for which `fcn` takes array parameters.

Example:
```
fcn(x) = x[1]^2 + (x[2]-1)^2
m = Minuit(fcn, [1, 0]; name = ["a", "b"], error = 0.1*ones(2), fix_a = true, limit_b = (0, 50) )
```
where the parameters are collected in an array `par` which is the argument of `fcn(par)`,
and `typeof(m) = ArrayFit`. In this case, one can use external code to compute the gradient as
`gradfun(par) = gradient(fcn, par)` (the exported `gradient` function is from `ForwardDiff`),
and include `grad = gradfun` as a keyword argument.

If `fcn` is defined as `fcn(a, b)`, then the starting values need to be set as
`Minuit(fcn, a = 1, b = 0)`.

`migrad, minos, hesse, matrix`:
wrappers of `iminuit.Minuit.migrad`, `iminuit.Minuit.minos`, `iminuit.Minuit.hesse`, `iminuit.Minuit.matrix`.



## Some useful functions

```
model_fit(model::Function, data::Data, start_values; kws...)
@model_fit model data start_values kws...
```
the returning type is `ArrayFit`, which can be passed
to `migrad`, `minos` etc. to perform the fit and error analysis.
 `model` is the function to be fitted to `data`; it should be of the form
`model(x, params)` with `params` given either as an array or a tuple.

`func_argnames(f::Function)`:  Extracting the argument names of a function as an array.

```
Data(x::T, y::T, err::T) where {T<:Vector{Real}}
Data(df::DataFrame)
```
Fields: `x, y, err, ndata`. This defines a type for data with three columns:` x, y, err`;
`ndata` is the number of data rows.

```
chisq(dist::Function, data, par; fitrange = ())
```
defines the χ² function: `fun` the function to be fitted to the data given by `data`.
The parameters are collected into `par`, given as an array or a tuple.
* `dist` should be defined as `dist(x, par)`.
* `data` can be either given as the `Data` type, or of the form `(xdata, ydata [, err])`.
If no `err` is given explicitly, the errors are assumed to be 1 for all data points.
* `fitrange`: default to the whole data set; may given as, e.g., `2:10`,
which means only fitting to the 2nd to the 10th data points.


```
@plt_data(data, kws...)
@plt_data!(data, kws...)
@plt_best(dist, fit, data, kws...)
```
convenient macros for plotting the data, and a comparison of the best-fit curve
with the data (needs `using Plots`). `kwds...` in `@plt_data` adjusts settings of the `scatter` plot, and that in `@plt_best` adjusts settings of the best-fit curve `plot`.


### Additional functions for error analysis

The following functions could be useful if `MINOS` fails to produce parameter uncertainties.

```
get_contours(fit::AbstractFit, χsq, parameters_combination::Vector{Int}; npts::Int=20, limits=true, sigma = 1.0)
```
For a given fit `fit` and the χ² function `χsq`, gives a list of parameter arrays,
with each array corresponding to a set of parameters obtained from calculating the `MINOS` 1σ contour
(try to find `npts` points in the contour) for the two parameters in `parameters_combination`.

`parameters_combination` is an `Int` array of the numbers of that two parameters, e.g. it is `[1, 2]` for the first
two parameters and `[2, 3]` or the second and third parameters.

If `limits` is `true`, then fix one parameter to its bounds from `MINOS` of the best fit and get the values for the other parameters using `MIGRAD`;
this runs over all parameters.


```
get_contours_all(fit::AbstractFit, χsq; npts=20, limits=true, sigma = 1.0)
```
For a given fit `fit` and the χ² function `χsq`, gives a list of parameters sets which are at the edge of
1σ `MINOS` contours for all parameters combinations. The case of `limits` being `true` runs only once.


`contour_df(fit::AbstractFit, χsq; npts=20, limits=true, sigma = 1.0)` gives such parameters in the form of a `DataFrame`.

```
get_contours_given_parameter(fit::AbstractFit, χsq, para::T, range) where {T <: Union{Symbol, String}}
contour_df_given_parameter(fit::AbstractFit, χsq, para::T, range; limits = true) where {T <: Union{Symbol, String}}
```
give parameter sets in 1σ for a given parameter constrained in a range (the latter returns a `DataFrame`).

For using array parameters, if no user-defined names have been given to the parameters,
`paras` should be given such that `"x0"` or `:x0` for the 1st parameter, `"x1"` or `:x1` for the 2nd parameter (default names in `iminuit`), ...


```
get_contours_samples(fit::AbstractFit, χsq, paras, ranges; nsamples = 100, MNbounds = true)
contour_df_samples(fit::AbstractFit, χsq, paras, ranges; nsamples = 100, MNbounds=true)
```
gives 1σ parameter sets as an `Array` (the latter returns a `DataFrame`) for given parameters constrained in `ranges`:
* if `paras` is a single parameter, then take equally spaced `nsamples` in `ranges` given in the form of `(min, max)`;
* if `paras` contain more (≧2) parameters, then `paras` should be of the form `(:para1, :para2)`, `ranges` should be of the form `((min1, max1), (min2, max2))`,
and values for the parameters given in `paras` are randomly sampled in the given `ranges`;
* if `MNbounds` is true, then constrain the parameters in the range provided by `MINOS` no matter whether that is valid or not (to be improved by checking the validity)
* if `igrad` is true, then use `ForwardDiff.gradient` to compute the gradient.
* For using array parameters, if no user-defined names have been given to the parameters,
`paras` should be given such that `"x0"` or `:x0` for the 1st parameter, `"x1"` or `:x1` for the 2nd parameter, ...
