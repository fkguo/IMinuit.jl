using PyCall
minuit = pyimport(:iminuit);

"""
    method_argnames(m::Method)

    Extracting the argument names of a method as an array.
    Modified from [https://github.com/JuliaLang/julia/blob/master/base/methodshow.jl]() (`Vector{Any}`` changed to `Vector{Symbol}`)
"""
function method_argnames(m::Method)
    argnames = ccall(:jl_uncompress_argnames, Vector{Symbol}, (Any,), m.slot_syms)
    isempty(argnames) && return argnames
    return argnames[1:m.nargs]
end

"""
    func_argnames(m::Method)

    Extracting the argument names of a function as an array.
"""
function func_argnames(f::Function)
    ms = collect(methods(f))
    return method_argnames(last(ms))[2:end]
end

# this defines only with kwds
"""
    Minuit(fcn; kwds...)
    Minuit(fcn, start; kwds...)

    Wrapper of the `iminuit` function `Minuit`.
    `fcn` is the function to be minized.
    `start`: an array/tuple of the starting values of the parameters.
    `kwds` is the list of keywrod arguments of `Minuit`. For more information, refer to the `iminuit` manual.

    Example:

    `Minuit(fcn, [1., 1, 0]; name = ["a", "b", "c"], error = 0.1*ones(3), fix_a = true, limit_b = (0, 50) )`
    where the parameters are collected in an array `par` which is the argument of `fcn(par)`. In this case,
    one can use external code (e.g., `using ForwardDiff: gradient`) to compute the gradient as `gradfun(par) = gradient(fcn, par)`, and include `grad = gradfun` as a keyword argument.


    `fcn` can also be defined as `fcn(a, b, c)`.

"""
function Minuit(fcn; kwds...)::PyObject
    forced_parameters = Tuple(func_argnames(fcn)) # get the argument lists of fcn
    return minuit.Minuit(fcn; forced_parameters= forced_parameters, pedantic = false, kwds...)
end

Minuit(fcn, start; kwds...)::PyObject =  minuit.Minuit.from_array_func(fcn, start; pedantic = false, kwds...)


function migrad(f::PyObject; ncall = 1000, resume = true, nsplit = 1, precision = nothing)
    return pycall(f.migrad, PyObject, ncall, resume, nsplit, precision)
end

migrad(f::PyObject) = pycall(f.migrad, PyObject)

hesse(f::PyObject; maxcall = 0) = pycall(f.hesse, PyObject, maxcall)
hesse(f::PyObject) = pycall(f.hesse, PyObject)

function minos(f::PyObject; var = nothing, sigma = 1, maxcall = 0)
    return pycall(f.minos, PyObject, var, sigma, maxcall)
end

minos(f::PyObject) = pycall(f.minos, PyObject)
