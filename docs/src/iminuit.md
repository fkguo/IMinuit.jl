# Functions in `iminuit`

## `Minuit`
```
Minuit(fcn; kwds...)::Fit
Minuit(fcn, start; kwds...)::ArrayFit
```

 Wrapper of the `iminuit` function `Minuit`.
 `fcn` is the function to be optimized.
 `start`: an array/tuple of the starting values of the parameters.
 `kwds` is the list of keywrod arguments of `Minuit`. For more information, refer to the [`iminuit` manual](https://iminuit.readthedocs.io/en/stable/about.html).
 The `Fit` one, for which `fcn` takes individual parameters as variables,  is generally
 faster than the `ArrayFit` one, for which `fcn` takes array parameters.


## `migrad`, `hesse`, `minos`

```
migrad, minos, hesse, matrix
```
wrappers of `iminuit.Minuit.migrad`, `iminuit.Minuit.minos`, `iminuit.Minuit.hesse`, `iminuit.Minuit.matrix`.


```
migrad(f::AbstractFit; ncall = 1000, resume = true, nsplit = 1, precision = nothing)
```
minimizing using `MIGRAD`. Further information can be found by
```@example
using IMinuit # hide
mMinuit = iminuit.Minuit # hide
for f in [:migrad, :minos, :hesse, :matrix, :args] # hide
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


```
hesse(f::AbstractFit; maxcall = 0)
```
run `HESSE` to compute parabolic errors.
```@example
using IMinuit # hide
mMinuit = iminuit.Minuit # hide
for f in [:migrad, :minos, :hesse, :matrix, :args] # hide
    sf = string(f) # hide
    @eval @doc LazyHelp(mMinuit, $sf)  function $f(ars...; kws...) #function $f(args...; kws...) # hide
       if !hasproperty(mMinuit, $sf) # hide
            error("iminuit ", version, " does not have iminuit.Minuit", $sf) # hide
        end # hide
        return pycall(mMinuit.$sf, PyAny, ars...; kws...) # hide
    end # hide
end # hide
@doc hesse
```


```
minos(f::AbstractFit; var = nothing, sigma = 1, maxcall = 0)
```
run `MINOS` to compute asymmetric errors taking into account correlation and nonlinearity
```@example
using IMinuit # hide
mMinuit = iminuit.Minuit # hide
for f in [:migrad, :minos, :hesse, :matrix, :args] # hide
    sf = string(f) # hide
    @eval @doc LazyHelp(mMinuit, $sf)  function $f(ars...; kws...) #function $f(args...; kws...) # hide
       if !hasproperty(mMinuit, $sf) # hide
            error("iminuit ", version, " does not have iminuit.Minuit", $sf) # hide
        end # hide
        return pycall(mMinuit.$sf, PyAny, ars...; kws...) # hide
    end # hide
end # hide
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



##  Some other useful `iminuit` functions

For a fit `f`, 
```
f.contour(par1, par2; kws)
f.mncontour(par1, par2; kws)
f.mnprofile(par1; kws)
```
```@example
using IMinuit # hide
mMinuit = iminuit.Minuit # hide
for f in [:contour] # hide
    sf = string(f) # hide
    @eval @doc LazyHelp(mMinuit, $sf)  function $f(ars...; kws...) #function $f(args...; kws...) # hide
       if !hasproperty(mMinuit, $sf) # hide
            error("iminuit ", version, " does not have iminuit.Minuit", $sf) # hide
        end # hide
        return pycall(mMinuit.$sf, PyAny, ars...; kws...) # hide
    end # hide
end # hide
@doc contour
```

```@example
using IMinuit # hide
mMinuit = iminuit.Minuit # hide
for f in [:mncontour] # hide
    sf = string(f) # hide
    @eval @doc LazyHelp(mMinuit, $sf)  function $f(ars...; kws...) #function $f(args...; kws...) # hide
       if !hasproperty(mMinuit, $sf) # hide
            error("iminuit ", version, " does not have iminuit.Minuit", $sf) # hide
        end # hide
        return pycall(mMinuit.$sf, PyAny, ars...; kws...) # hide
    end # hide
end # hide
@doc mncontour
```

```@example
using IMinuit # hide
mMinuit = iminuit.Minuit # hide
for f in [:mnprofile] # hide
    sf = string(f) # hide
    @eval @doc LazyHelp(mMinuit, $sf)  function $f(ars...; kws...) #function $f(args...; kws...) # hide
       if !hasproperty(mMinuit, $sf) # hide
            error("iminuit ", version, " does not have iminuit.Minuit", $sf) # hide
        end # hide
        return pycall(mMinuit.$sf, PyAny, ars...; kws...) # hide
    end # hide
end # hide
@doc mnprofile
```

Other `iminuit` functions can also be used, like
`f.draw_contour(par1, par2; kws)` and `f.draw_mncontour(par1, par2; kws)`.