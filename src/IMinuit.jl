__precompile__() # this module is safe to precompile
module IMinuit

using PyCall: PyObject, pycall, PyNULL, PyAny, PyVector, pyimport_conda, pyimport, pytype_mapping
import PyCall: PyObject, pycall
# import PyCall: hasproperty # Base.hasproperty in Julia 1.2
import Base: convert, ==, isequal, hash, hasproperty,  haskey

using ForwardDiff: gradient

export Minuit, migrad, minos, hesse, matrix, iminuit, args, model_fit, @model_fit
export AbstractFit, Fit, ArrayFit, func_argnames, Data, chisq, @plt_data, @plt_data!, @plt_best, @plt_best!
export gradient, LazyHelp
export get_contours, get_contours_all, contour_df, get_contours_given_parameter
export contour_df_given_parameter, get_contours_samples, contour_df_samples

# copied from PyPlot.jl
# that lazily looks up help from a PyObject via zero or more keys.
# This saves us time when loading iminuit, since we don't have
# to load up all of the documentation strings right away.
struct LazyHelp
    o # a PyObject or similar object supporting getindex with a __doc__ property
    keys::Tuple{Vararg{String}}
    LazyHelp(o) = new(o, ())
    LazyHelp(o, k::AbstractString) = new(o, (k,))
    LazyHelp(o, k1::AbstractString, k2::AbstractString) = new(o, (k1,k2))
    LazyHelp(o, k::Tuple{Vararg{AbstractString}}) = new(o, k)
end
function show(io::IO, ::MIME"text/plain", h::LazyHelp)
    o = h.o
    for k in h.keys
        o = getproperty(o, k) # o[k]
    end
    if hasproperty(o, "__doc__")
        print(io, convert(AbstractString, o."__doc__"))
    else
        print(io, "no Python docstring found for ", h.k)
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


"""
    method_argnames(m::Method)

    Extracting the argument names of a method as an array.
    Modified from [`methodshow.jl`](https://github.com/JuliaLang/julia/blob/master/base/methodshow.jl) (`Vector{Any}`` changed to `Vector{Symbol}`)
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

###########################################################################


# Wrappers of the iminuit functions
@doc raw"""
    Minuit(fcn; kwds...)
    Minuit(fcn, start; kwds...)

Wrapper of the `iminuit` function `Minuit`.
* `fcn` is the function to be optimized.
* `start`: an array/tuple of the starting values of the parameters.
* `kwds` is the list of keywrod arguments of `Minuit`. For more information, refer to the `iminuit` manual.

Example:
```
fit = Minuit(fcn, [1, 0]; name = ["a", "b"], error = 0.1*ones(2), 
            fix_a = true, limit_b = (0, 50) )
migrad(fit)
```
where the parameters are collected in an array `par` which is the argument of `fcn(par)`. In this case,
one can use external code (e.g., `using ForwardDiff: gradient`) to compute the gradient as `gradfun(par) = gradient(fcn, par)`, and include `grad = gradfun` as a keyword argument.


If `fcn` is defined as `fcn(a, b)`, then the starting values need to be
set as `Minuit(fcn, a = 1, b = 0)`.

From `iminuit`:

`Minuit(fcn, throw_nan=False, pedantic=True, forced_parameters=None, print_level=0, 
        errordef=None, grad=None, use_array_call=False, **kwds)`

"""
function Minuit(fcn; kwds...)::Fit
    forced_parameters = Tuple(func_argnames(fcn)) # get the argument lists of fcn
    return iminuit.Minuit(fcn; forced_parameters= forced_parameters, pedantic = false, kwds...)
end

Minuit(fcn, start::AbstractVector; kwds...)::ArrayFit =  iminuit.Minuit.from_array_func(fcn, start; pedantic = false, kwds...)
Minuit(fcn, start::Tuple; kwds...)::ArrayFit =  iminuit.Minuit.from_array_func(fcn, start; pedantic = false, kwds...)

function migrad(f::AbstractFit; ncall = 1000, resume = true, nsplit = 1, precision = nothing)
    return pycall(f.migrad, PyObject, ncall, resume, nsplit, precision)
end


hesse(f::AbstractFit; maxcall = 0) = pycall(f.hesse, PyObject, maxcall)


function minos(f::AbstractFit; var = nothing, sigma = 1, maxcall = 0)
    return pycall(f.minos, PyObject, var, sigma, maxcall)
end


matrix(f::AbstractFit; kws...) = pycall(PyObject(f).matrix, PyObject, kws...)

function args(o::AbstractFit)::Vector{Float64}
    _a::PyObject = o.args;  _n::Int = length(_a)
    _res = zeros(_n)
    @views for i = 1:_n
        _res[i] = get(_a, i-1)
    end
    return _res
end

# for f in [:migrad, :minos, :hesse, :matrix, :args, :contour, :mncontour, :mnprofile]
#     sf = string(f)
#     @eval @doc LazyHelp(mMinuit, $sf)  function $f(ars...; kws...) #function $f(args...; kws...)
#         if !hasproperty(mMinuit, $sf)
#             error("iminuit ", version, " does not have iminuit.Minuit", $sf)
#         end
#         return pycall(mMinuit.$sf, PyAny, ars...; kws...)
#     end
# end

#########################################################################

include("Data.jl")
include("contour.jl")

#########################################################################

"""
    model_fit(model::Function, data::Data, start_values; kws...)

convenient wrapper, the returning stype is `ArrayFit`, which can be passed
to `migrad`, `minos` etc.
* `model` is the function to be fitted to `data`; it should be of the form
`model(x, params)` with `params` given either as an array or a tuple.
"""
function model_fit(model::Function, data::Data, start_values; kws...)
    # _model(x, par) = length(func_argnames(model)) > 2 ? model(x, par...) : model(x, par)
    _chisq(par) = chisq(model, data, par)
    _fit = Minuit(_chisq, start_values; kws...)
    return _fit
end

"""
    @model_fit model data start_values kws...

convenient wrapper, the returning stype is `ArrayFit`, which can be passed
to `migrad`, `minos` etc.
* `model` is the function to be fitted to `data`; it should be of the form
`model(x, params)` with `params` given either as an array or a tuple.
"""
macro model_fit(model, data, start_values, kws...)
    _expr = quote
        _chisq(par) = chisq($model, $data, par)
        _fit = isempty($kws) ? Minuit(_chisq, $start_values) : Minuit(_chisq, $start_values; $(kws...))
    end
    esc(_expr)
end

end
