__precompile__() # this module is safe to precompile
module IMinuit

using PyCall: PyObject, pycall, PyNULL, PyAny, PyVector, pyimport_conda, pyimport, pytype_mapping, set!
using Distributed
using LinearAlgebra: norm
# import PyCall: hasproperty # Base.hasproperty in Julia 1.2
import Base: convert, ==, isequal, hash, hasproperty, haskey
import PyCall.pycall

using ForwardDiff: gradient

export Minuit, migrad, minos, hesse, iminuit, args, model_fit, @model_fit
export AbstractFit, Fit, ArrayFit, func_argnames, Data, chisq, @plt_data, @plt_data!, @plt_best, @plt_best!
export gradient, LazyHelp, contour, mncontour, mnprofile, draw_mncontour, profile,
    draw_contour, draw_profile
export matrix
export get_contours, get_contours_all, contour_df, get_contours_given_parameter
export contour_df_given_parameter, get_contours_samples, contour_df_samples
export reset
export chi2, poisson_chi2, multinominal_chi2
export scipy, simplex, scan, draw_mnmatrix

# copied from PyPlot.jl
# that lazily looks up help from a PyObject via zero or more keys.
# This saves us time when loading iminuit, since we don't have
# to load up all of the documentation strings right away.
struct LazyHelp
    o # a PyObject or similar object supporting getindex with a __doc__ property
    keys::Tuple{Vararg{String}}
    LazyHelp(o) = new(o, ())
    LazyHelp(o, k::AbstractString) = new(o, (k,))
    LazyHelp(o, k1::AbstractString, k2::AbstractString) = new(o, (k1, k2))
    LazyHelp(o, k::Tuple{Vararg{AbstractString}}) = new(o, k)
end
# function show(io::IO, ::MIME"text/plain", h::LazyHelp)
# 	o = h.o
#   for k in h.keys
#   	o = getproperty(o, k)
#   end
#   if hasproperty(o, "__doc__")
#   	print(io, "Docstring pulled from the Python `iminuit`:\n\n")
#     print(io, convert(AbstractString, o.))
#   end
# end
function show(io::IO, ::MIME"text/plain", h::LazyHelp)
    o = h.o
    for k in h.keys
        o = getproperty(o, k)
    end
    if hasproperty(o, "__doc__")
        print(io, "Docstring pulled from the Python `iminuit`:\n\n")
        print(io, convert(AbstractString, o."__doc__"))
    else
        print(io, "no Python docstring found for ", h.keys)
    end
end

Base.show(io::IO, h::LazyHelp) = show(io, "text/plain", h)
function Base.Docs.catdoc(hs::LazyHelp...)
    Base.Docs.Text() do io
        for h in hs
            show(io, MIME"text/plain"(), h)
        end
    end
end


@doc raw"""
method_argnames(m::Method)

Extracting the argument names of a method as an array.
Modified from [`methodshow.jl`](https://github.com/JuliaLang/julia/blob/master/base/methodshow.jl) 
(`Vector{Any}` changed to `Vector{Symbol}`).
"""
function method_argnames(m::Method)
    argnames = ccall(:jl_uncompress_argnames, Vector{Symbol}, (Any,), m.slot_syms)
    isempty(argnames) && return argnames
    return argnames[1:m.nargs]
end

"""
func_argnames(f::Function)

Extracting the argument names of a function as an array.
"""
function func_argnames(f::Function)
    ms = collect(methods(f))
    return method_argnames(last(ms))[2:end]
end


###########################################################################

include("init.jl")
include("FitStructs.jl")
include("preprocess.jl")

###########################################################################


# Wrappers of the iminuit functions
@doc raw"""
Minuit(fcn; kwds...)
Minuit(fcn, start; kwds...)
Minuit(fcn, m::AbatractFit; kwds...)

> In the v2.0 or later, many keywords are removed and replaced by new usage, in this package you can still use the old usage like passing `limit_`, `fix_` and the like to the `Minuit` function, but the code to implement these usage is over complex, if you have better implementation you can pull a request, or you can learn the new usage if you want to keep your code robust and reusable

Wrapper of the `iminuit` function `Minuit`.
* `fcn` is the function to be optimized.
* `start`: an array/tuple of the starting values of the parameters.
* `kwds` is the list of keyword arguments of `Minuit`. For more information, refer to the `iminuit` manual.
* `m`: a fit that was previously defined; its parameters at the latest stage can be passed to the new fit.

Example:
```
fit = Minuit(fcn, [1, 0]; name = ["a", "b"], error = 0.1*ones(2), 
fix_a = true, limit_b = (0, 50) )
migrad(fit)
```
where the parameters are collected in an array `par` which is the argument of `fcn(par)`. In this case,
one can use external code (e.g., `using ForwardDiff: gradient`) to compute the gradient as 
`gradfun(par) = gradient(fcn, par)`, and include `grad = gradfun` as a keyword argument.


If `fcn` is defined as `fcn(a, b)`, then the starting values need to be
set as `Minuit(fcn, a = 1, b = 0)`.

From `iminuit`:

`Minuit(fcn: Callable, *args: Union[float, Sequence[float]], grad: Callable = None, name: Collection[str] = None, **kwds: float)`

"""
function Minuit(fcn; kwds...)::Fit
    removed = ["errordef", "throw_nan", "print_level", "use_array_call"]
    new_kwds, stored_kwds, fitarg_dict = preprocess(fcn; kwds...)
    # println(new_kwds)
    m = iminuit.Minuit(fcn; new_kwds...)

    if !all(fitarg_dict["error"] .== 0.0)
        m.errors = fitarg_dict["error"]
    end
    m.limits = fitarg_dict["limit"]
    m.fixed = fitarg_dict["fix"]

    for sym in removed
        if haskey(stored_kwds, sym)
            m[Symbol(sym)] = stored_kwds[sym]
        end
    end
    return m
end

function Minuit(fcn, start::AbstractVector; kwds...)::ArrayFit
    removed = ["errordef", "throw_nan", "print_level", "use_array_call"]
    if isempty(kwds)
        m = iminuit.Minuit(fcn, start)
    else
        if !haskey(kwds, :name)
            name = collect(Symbol(:x, i - 1) for i in 1:length(start))
            # println(name)
            kwds = [:name => name, kwds...]
        end
        # println(kwds)
        new_kwds, stored_kwds, fitarg_dict = preprocess(fcn; kwds...)
        m = iminuit.Minuit(fcn, start; new_kwds...)

        if !all(fitarg_dict["error"] .== 0.0)
            m.errors = fitarg_dict["error"]
        end
        m.limits = fitarg_dict["limit"]
        m.fixed = fitarg_dict["fix"]

        for sym in removed
            if haskey(stored_kwds, sym)
                m[Symbol(sym)] = stored_kwds[sym]
            end
        end
    end
    return m
end

function Minuit(fcn, start::Tuple; kwds...)::ArrayFit
    removed = ["errordef", "throw_nan", "print_level", "use_array_call"]
    if isempty(kwds)
        m = iminuit.Minuit(fcn, start)
    else
        new_kwds, stored_kwds, fitarg_dict = preprocess(fcn; kwds...)
        m = iminuit.Minuit(fcn, start; new_kwds...)

        if !all(fitarg_dict["error"] .== 0.0)
            m.errors = fitarg_dict["error"]
        end
        m.limits = fitarg_dict["limit"]
        m.fixed = fitarg_dict["fix"]

        for sym in removed
            if haskey(stored_kwds, sym)
                m[Symbol(sym)] = stored_kwds[sym]
            end
        end
    end
    return m
end

function Minuit(fcn, m::AbstractFit; kwds...)
    removed = ["errordef", "throw_nan", "print_level", "use_array_call"]
    new_kwds, stored_kwds, ini_value, fitarg_dict = preprocess(fcn, m; kwds...)
    _m::ArrayFit = Minuit(fcn, ini_value; new_kwds...)

    if !all(fitarg_dict["error"] .== 0.0)
        _m.errors = fitarg_dict["error"]
    end
    _m.limits = fitarg_dict["limit"]
    _m.fixed = fitarg_dict["fix"]

    for sym in removed
        if haskey(stored_kwds, sym)
            m[Symbol(sym)] = stored_kwds[sym]
        end
    end
    return _m
end

reset(f::AbstractFit) = pycall(f.reset, PyObject)

set_precision(f::AbstractFit, new_precision) = set!(PyObject(f), precision, new_precision)

# nsplit = 1, nsplit option removed since iminuit v1.5.0
# resume was removed, use reset() instead
# precision was also removed, use set_precision() instead
function migrad(f::AbstractFit; ncall=nothing, resume=true, precision=nothing)
    if !resume
        reset(f)
    end
    if precision !== nothing
        f.precision = precision
    end
    return pycall(f.migrad, PyObject, ncall)
end



hesse(f::AbstractFit; maxcall=0) = pycall(f.hesse, PyObject, maxcall)

function minos(f::AbstractFit, var=nothing, sigma=1, maxcall=0)
    if sigma >= 1
        if var === nothing
            return pycall(f.minos, PyObject, cl=sigma, ncall=maxcall)
        else
            return pycall(f.minos, PyObject, (var), sigma, maxcall)
        end
    else
        error("Sigma less than 1, not supported yet!")
    end
end

# matrix(f::AbstractFit; kws...) = pycall(PyObject(f).matrix, PyObject, kws...)
function args(o::AbstractFit)::Vector{Float64}
    _a::PyObject = o.values
    _n::Int = length(_a)
    _res = zeros(_n)
    @views for i = 1:_n
        _res[i] = get(_a, i - 1)
    end
    return _res
end

mncontour(f::AbstractFit, par1, par2; numpoints=100, sigma=nothing, kws...) = f.mncontour(par1, par2; cl = sigma, size = numpoints, kws...)
draw_mncontour(f::AbstractFit, par1, par2; numpoints=100, nsigma=nothing, kws...) = f.draw_mncontour(par1, par2; cl = nsigma, size = numpoints, kws...)

for fun in [:scan, :simplex]
    :(($fun)(f::AbstractFit, par1; kws...) = f.$fun(par1; kws...)) |> eval
end
#fix for incorrect parameters
for fun in [:mnprofile, :draw_mnprofile]
    :(($fun)(f::AbstractFit, par; bins = 30, kws...) = begin
        f.$fun(par; size = bins, kws...)
    end) |> eval
end

for fun in [:profile, :draw_profile]
    :(($fun)(f::AbstractFit, par; bins = 100, kws...) = begin
        # println("usage of bin is deprecated, please use the arguments size, for more info, go check the iminuit python document")
        f.$fun(par; size = bins, kws...)
    end) |> eval
end

for fun in [:contour, :draw_contour]
    :(($fun)(f::AbstractFit, par1, par2; bins = 50, kws...) = begin
        f.$fun(par1, par2; size = bins, kws...)
    end) |> eval
end

for fun in [:scipy, :draw_mnmatrix]
    :(($fun)(f::AbstractFit; kws...) = begin
        f.$fun(kws...)
    end) |> eval
end


for f in [:migrad, :minos, :hesse, :contour, :mncontour, :profile,
    :mnprofile, :draw_mncontour, :draw_contour, :draw_profile, :draw_mnprofile, :scan, :simplex, :scipy, :draw_mnmatrix]
    sf = string(f)
    @eval @doc LazyHelp(mMinuit, $sf) function $f(ars...; kws...)
        if !hasproperty(mMinuit, $sf)
            error("iminuit ", iminuit_version, " does not have iminuit.Minuit", $sf)
        end
        return pycall(mMinuit.$sf, PyAny, ars...; kws...)
    end
end

for f in [:chi2, :multinominal_chi2, :poisson_chi2]
    sf = string(f)
    @eval @doc LazyHelp(cost, $sf) function $f(ars...; kws...)
        if !hasproperty(cost, $sf)
            error("iminuit ", iminuit_version, " does not have iminuit.cost", $sf)
        end
        return pycall(cost.$sf, PyAny, ars...; kws...)
    end
end

#########################################################################

include("Data.jl")
include("contour.jl")

#########################################################################


"""
model_fit(model::Function, data::Data, start_values; kws...)

convenient wrapper for fitting a `model` to `data`; the returning stype is `ArrayFit`, which can be passed
to `migrad`, `minos` etc.
* `model` is the function to be fitted to `data`; it should be of the form `model(x, params)` with `params` given either as an array or a tuple.
* `start_values` can be either the initial values of parameters (`<: AbstractArray` or `<:Tuple`) or a previous fit of type `AbstractFit`.
"""
function model_fit(model::Function, data::Data, start_values; kws...)
    # _model(x, par) = length(func_argnames(model)) > 2 ? model(x, par...) : model(x, par)
    _chisq(par) = chisq(model, data, par)
    _fit = Minuit(_chisq, start_values; kws...)
    return _fit
end
function model_fit(model::Function, data::Data, fit::AbstractFit; kws...)
    _chisq(par) = chisq(model, data, par)
    _fit = Minuit(_chisq, fit; kws...)
    return _fit
end

"""
@model_fit model data start_values kws...

convenient wrapper for fitting a `model` to `data`; see [`model_fit`](@ref).
"""
macro model_fit(model, data, start_values, kws...)
    _expr = quote
        _chisq(par) = chisq($model, $data, par)
        _fit = isempty($kws) ? Minuit(_chisq, $start_values) : Minuit(_chisq, $start_values; $(kws...))
    end
    esc(_expr)
end

function matrix(f::AbstractFit; correlation=false, skip_fixed=true)
    if !f.valid
        f.migrad()
    end
    if skip_fixed == true
        if correlation == true
            return PyObject(f)."covariance".correlation()
        else
            return f.covariance
        end
    else
        error("the case that skip_fixed is false hasn't been implemented yet!")
    end
end


end
