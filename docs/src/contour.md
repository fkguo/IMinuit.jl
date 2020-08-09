# Functions for error analysis

There are functions in `iminuit` (`MINUIT`) to analyze parameter errors. Suppose we have an `AbstractFit` object `fit`, then one can use functions like `fit.draw_contour(:par1, :par2)`, `fit.draw_mncontour(:par1, :par2)` etc., see the [`iminit` documentation](https://iminuit.readthedocs.io/en/stable/reference.html).


The following extra functions could also be useful to get a set of parameters for error analysis if the model is highly nonlinear in the parameters.

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
1σ `MINOS` contours for all combinations of varying parameters. The case of `limits` being `true` runs only once.


`contour_df(fit::AbstractFit, χsq; npts=20, limits=true, sigma = 1.0)`
returns such parameters in the form of a `DataFrame`.

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
