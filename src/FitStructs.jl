abstract type AbstractFit end

# # modified from `Figure` in PyPlot
"""
    Fit <: AbstractFit

struct for fit with individual parameters.
"""
mutable struct Fit <: AbstractFit
    o::PyObject
end

"""
    ArrayFit <: AbstractFit

struct for fit with parameters collected in an `Array`.
"""
mutable struct ArrayFit <: AbstractFit
    o::PyObject
end

# modified from `Figure` in PyPlot
PyObject(f::T) where {T<:AbstractFit} = getfield(f, :o)
convert(::Type{T}, o::PyObject) where {T<:AbstractFit} = T(o)
==(f::T, g::T) where {T<:AbstractFit} = PyObject(f) == PyObject(g)
==(f::AbstractFit, g::PyObject) = PyObject(f) == g
==(f::PyObject, g::AbstractFit) = f == PyObject(g)
hash(f::AbstractFit) = hash(PyObject(f))
pycall(f::AbstractFit, args...; kws...) = pycall(PyObject(f), args...; kws...)
(f::AbstractFit)(args...; kws...) = pycall(PyObject(f), PyAny, args...; kws...)
Base.Docs.doc(f::AbstractFit) = Base.Docs.doc(PyObject(f))
#
Base.setproperty!(f::AbstractFit, s::Symbol, x) = setproperty!(PyObject(f), s, x)
Base.setproperty!(f::AbstractFit, s::AbstractString, x) = setproperty!(PyObject(f), s, x)
hasproperty(f::AbstractFit, s::Symbol) = hasproperty(PyObject(f), s)
hasproperty(f::AbstractFit, s::AbstractString) = hasproperty(PyObject(f), s)
Base.propertynames(f::AbstractFit) = propertynames(PyObject(f))
haskey(f::AbstractFit, x) = haskey(PyObject(f), x)

# let the matrix property to be given as an array
# matrix removed, now using Minuit.covariance.correlation(), might have some problem cause the return type is slightly different
function Base.getproperty(f::AbstractFit, s::Symbol)
    # use PyObject(f) instead of f to prevent StackOverflowError
    s === :matrix ? matrix(f) : getproperty(PyObject(f), s)
end
function Base.getproperty(f::AbstractFit, s::AbstractString)
    s === "matrix" ? matrix(f) : getproperty(PyObject(f), s)
end

