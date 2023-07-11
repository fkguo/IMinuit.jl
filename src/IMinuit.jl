__precompile__() # this module is safe to precompile
module IMinuit

using PyCall: PyObject, pycall, PyNULL, PyAny, PyVector, pyimport_conda, pyimport, pytype_mapping, set!
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
    print(io, convert(AbstractString, o.\"__doc__"))
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

###########################################################################


# Wrappers of the iminuit functions
@doc raw"""
Minuit(fcn; kwds...)
Minuit(fcn, start; kwds...)
Minuit(fcn, m::AbatractFit; kwds...)

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

`Minuit(fcn, throw_nan=False, pedantic=True, forced_parameters=None, print_level=0, 
errordef=None, grad=None, use_array_call=False, **kwds)`

"""
fitarg = ["error", "fix", "limit"]
removed = ["pedantic", "errordef", "throw_nan", "print_level", "use_array_call"]
removed_syms = [:pedantic, :errordef, :throw_nan, :print_leve, :use_array_call]
function preprocess(fcn; kwds...)
  args = nothing
  if haskey(kwds, :name)
  	args = kwds[:name]
  else
    args = func_argnames(fcn)
  end
  # println(args)
  arg_dict = Dict(x => i for (i, x) in enumerate(args))
  len = length(args)
  stored_kwds = Dict(str => nothing for str in removed)
  new_kwds = Vector()
  foo = Dict{String, Array{Union{Float64, Bool, Tuple{Float64, Float64}}}}()
  foo["error"] = Array{Float64}(undef, len)
  foo["limit"] = fill((-Inf64, Inf64), len)
  foo["fix"] = falses(len)
  for (k, v) in kwds
    k_str = String(k)
    if k_str in removed
    	stored_kwds[k_str] = v
      continue
    elseif k_str in fitarg
      foo[k_str] = v
      continue
    end
    udscore = findfirst('_', k_str)
    if udscore === nothing
      push!(new_kwds, (k, v))
    	continue
    end
    typ = k_str[1: udscore - 1]
    if typ in fitarg
      para = Symbol(k_str[udscore + 1:end])
      foo[typ][arg_dict[para]] = v
    end
  end
  return (new_kwds, foo)
end
function Minuit(fcn; kwds...)::Fit
#   forced_parameters = Tuple(func_argnames(fcn)) # get the argument lists of fcn
  # return iminuit.Minuit(fcn; forced_parameters=forced_parameters, kwds...)
  # args = func_argnames(fcn)
  # len = length(args)
  # arg_dict = Dict(x => i for (i, x) in enumerate(args))
  # stored_kwds = Dict(str => nothing for str in removed)
  # new_kwds = Vector()
  # foo = Dict{String, Array{Union{Float64, Bool, Tuple{Float64, Float64}}}}()
  # foo["error"] = Array{Float64}(undef, len)
  # foo["limit"] = fill((-Inf64, Inf64), len)
  # foo["fix"] = falses(len)
  # for (k, v) in kwds
  #   k_str = String(k)
  #   if k_str in removed
  #   	stored_kwds[k_str] = v
  #     continue
  #   end
  #   udscore = findfirst('_', k_str)
  #   if udscore === nothing
  #     push!(new_kwds, (k, v))
  #   	continue
  #   end
  #   typ = k_str[1: udscore - 1]
  #   if typ in fitarg
  #     para = Symbol(k_str[udscore + 1:end])
  #     foo[typ][arg_dict[para]] = v
  #   end
  # end
  new_kwds, foo = preprocess(fcn; kwds...)
  m = iminuit.Minuit(fcn; new_kwds...)
  m.errors = foo["error"]
  m.limits = foo["limit"]
  m.fixed = foo["fix"]
  # @remo removed_syms
  # for (kwd, v) in stored_kwds
  #   if v !== nothing
  #     :(m.$kwd = v) |> eval
  #   end
  # end
  return m
end

function Minuit(fcn, start::AbstractVector; kwds...)::ArrayFit
  new_kwds, foo = preprocess(fcn; kwds...)
  m = iminuit.Minuit(fcn, start; new_kwds...)
  m.errors = foo["error"]
  m.limits = foo["limit"]
  m.fixed = foo["fix"]
  return m
end
function Minuit(fcn, start::Tuple; kwds...)::ArrayFit
  new_kwds, foo = preprocess(fcn; kwds...)
  m = iminuit.Minuit(fcn, start; new_kwds...)
  m.errors = foo["error"]
  m.limits = foo["limit"]
  m.fixed = foo["fix"]
  return m
end
# Minuit(fcn; kwds...)::ArrayFit = iminuit.Minuit(fcn; kwds...)

# Minuit(fcn, m::ArrayFit; kwds...) = Minuit(fcn, args(m); (Symbol(k) => v for (k, v) in m.fitarg)..., kwds...)
function Minuit(fcn, m::ArrayFit; kwds...)
  # attributes = [:errors, :fixed, :limits]
  _m = Minuit(fcn; kwds...)
  _m.errors = m.errors
  _m.fixed = m.fixed
  _m.limits = m.limits
  # for attri in attributes
  # 	:(_m.$attri = __m.$attri) |> eval
  # end
  return _m
end
# Minuit(fcn, m::Fit; kwds...) = Minuit(fcn; (Symbol(k) => v for (k, v) in m.fitarg)..., kwds...)
function Minuit(fcn, m::Fit; kwds...)
  # attributes = [:errors, :fixed, :limits]
  _m = Minuit(fcn; kwds...)
  _m.errors = m.errors
  _m.fixed = m.fixed
  _m.limits = m.limits
  # for attri in attributes
  # 	:(_m.$attri = m.$attri) |> eval
  # end
  return _m
end

reset(f::AbstractFit) = pycall(f.reset, PyObject)

set_precision(f::AbstractFit, new_precision) = set!(PyObject(f), precision, new_precision)

# nsplit = 1, nsplit option removed since iminuit v1.5.0
# resume was removed, use reset() instead
# precision was also removed, use set_precision() instead
function migrad(f::AbstractFit; ncall=1000, resume=true, precision=nothing)
  if !resume
    reset(f)
  end
  if precision !== nothing
    set_precision(precision)
  end
  return pycall(f.migrad, PyObject, ncall)
end



hesse(f::AbstractFit; maxcall=0) = pycall(f.hesse, PyObject, maxcall)

function minos(f::AbstractFit)
  return pycall(f.minos, PyObject)
end

# new usage, now minos() can take confidence level as parameter, which should be less than 1
# function minos(f::AbstractFit, parameter = nothing, cl=nothing, maxcall=0)
#   return pycall(f.minos, PyObject, parameter, cl, maxcall)
# end

# legacy usage, pass the number of sigma into the function, which is exactly the same as the new usage, cause when parameter cl is larger than 1 it is interpreted as the number of sigma
# while the scenario when the number is less than 1 is not implemented, hope u guys do not need it LOL.
function minos(f::AbstractFit, var = nothing, sigma = 1, maxcall = 0)
  if sigma >= 1
    if var === nothing
      return pycall(f.minos, PyObject, cl = sigma, ncall = maxcall)
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

for fun in [:contour, :mncontour, :draw_contour, :draw_mncontour]
  :(($fun)(f::AbstractFit, par1, par2; kws...) = f.$fun(par1, par2; kws...)) |> eval
end
#fix for incorrect parameters
for fun in [:mncontour, :draw_mncontour]
  :(($fun)(f::AbstractFit, par1, par2; numpoints = 100, sigma = 1, kws...) = begin println("usage of numpoints or sigma is deprecated, please use the arguments size and cl, for more info, go check the iminuit python document"); f.$fun(par1, par2; cl = sigma, size = numpoints, kws...) end) |> eval
end

for fun in [:profile, :draw_profile, :mnprofile, :draw_mnprofile]
  :(($fun)(f::AbstractFit, par1; kws...) = f.$fun(par1; kws...)) |> eval
end
#fix for incorrect parameters
for fun in [:profile, :draw_profile, :mnprofile, :draw_mnprofile]
  :(($fun)(f::AbstractFit, par1; numpoints = 100, sigma = 1, kws...) = begin println("usage of numpoints or sigma is deprecated, please use the arguments size and cl, for more info, go check the iminuit python document"); f.$fun(par1; cl = sigma, size = numpoints, kws...) end) |> eval
end


for f in [:migrad, :minos, :hesse, :contour, :mncontour, :profile,
          :mnprofile, :draw_mncontour, :draw_contour, :draw_profile, :draw_mnprofile]
  sf = string(f)
  @eval @doc LazyHelp(mMinuit, $sf) function $f(ars...; kws...)
    if !hasproperty(mMinuit, $sf)
      error("iminuit ", iminuit_version, " does not have iminuit.Minuit", $sf)
    end
    return pycall(mMinuit.$sf, PyAny, ars...; kws...)
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

function matrix(f::AbstractFit; correlation = true, skip_fixed = true)
  if !f.valid
  	f.migrad()
  end
	if skip_fixed == true
    if correlation == true 
      return PyObject(f)."covariance".correlation()
    else
      error("still don't know how to implement the error matrix!")
    end
  else
    error("the case that skip_fixed is false hasn't been implemented yet!")
  end
end


end
