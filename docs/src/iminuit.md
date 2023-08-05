# Functions in `iminuit`

## `Minuit`
```
Minuit(fcn; kwds...)::Fit
Minuit(fcn, start; kwds...)::ArrayFit
Minuit(fcn, m::AbatractFit; kwds...)
```

 Wrapper of the `iminuit` function `Minuit`; for more information, refer to the [`iminuit` manual](https://iminuit.readthedocs.io/en/stable/about.html).
 * `fcn` is the function to be optimized.
 * `start`: an array/tuple of the starting values of the parameters.
 * `kwds` is the list of keyword arguments of `Minuit`. 
 * For the `Fit` one, `fcn` takes individual parameters as variables; for the `ArrayFit` one, `fcn` takes array parameters.
 * `m`: a fit (either `Fit` or `ArrayFit`) that was previously defined; its parameters at the latest stage can be passed to the new fit.


## `migrad`, `hesse`, `minos`

```
migrad, minos, hesse, matrix
```
wrappers of `iminuit.Minuit.migrad`, `iminuit.Minuit.minos`, `iminuit.Minuit.hesse`, `iminuit.Minuit.matrix`.

```
migrad(f::AbstractFit; ncall = nothing, resume = true, precision = nothing)
```
minimizing using `MIGRAD`. Further information (docstring from `iminuit`) can be found by
```@example 1
using IMinuit # hide
mMinuit = iminuit.Minuit # hide
for f in [:migrad, :minos, :hesse, :contour, :mncontour, :profile, :mnprofile] # hide
    sf = string(f) # hide
    @eval @doc LazyHelp(mMinuit, $sf)  function $f(ars...; kws...) #function $f(args...; kws...) # hide
       if !hasproperty(mMinuit, $sf) # hide
            error("iminuit ", version, " does not have iminuit.Minuit", $sf) # hide
        end # hide
        return pycall(mMinuit.$sf, PyAny, ars...; kws...) # hide
    end # hide
end # hide
@doc migrad
```
Note that the python package doesn't provide parameters like `resume` and `precision`, for which are implemented in this julia package.

```
hesse(f::AbstractFit; maxcall = 0)
```
run `HESSE` to compute parabolic errors.
```@example 1
@doc hesse
```


```
minos(f::AbstractFit; var = nothing, sigma = 1, maxcall = 0)
```
run `MINOS` to compute asymmetric errors taking into account correlation and nonlinearity
```@example 1
@doc minos
```


For a fit `f`, 
```
args(f::AbstractFit)
```
return the parameter values from fit `f`.

```
matrix(f::AbstractFit; kws...) 
```
return error or correlation matrix (set the keyword `correlation = true`).

### Other algorithms that can be used to minimize
For a fit `f`, 
```
f.scan(ncall)
f.scipy(kws...)
f.simplex(ncall)
```
> `scipy` method need the corresponding python package to be installed in the enviroment that `PyCall.jl` was linked to
```@example 1
@doc scan
```
```@example 1
@doc scipy
```
```@example 1
@doc simplex
```

##  Some other useful `iminuit` functions

For a fit `f`, 
```
f.contour(par1, par2; kws)
f.mncontour(par1, par2; kws)
f.draw_contour(par1, par2; kws)
f.draw_mncontour(par1, par2; kws)
f.profile(par1; kws)
f.draw_profile(par1; kws)
```

They can also be written as `contour(f, par1, par2; kws)` etc. Some docstrings from `iminuit`:

```@example 1
@doc contour
```

```@example 1
@doc mncontour
```

```@example 1
@doc mnprofile
```

Other `iminuit` functions can also be used, like
`f.draw_contour(par1, par2; kws)` and `f.draw_mncontour(par1, par2; kws)`.