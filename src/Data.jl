using DataFrames: DataFrame

"""
    Data(x::T, y::T, err::T) where {T<:Vector{Real}}
    Data(df::DataFrame)

Fields: `x, y, err, ndata`

This defines a type for data with three columns:` x, y, err`; `ndata` is the number of data rows.

Only symmetric errors (of `y`) are supported.
"""
struct Data
    x::Vector{Real}
    y::Vector{Real}
    err::Vector{Real}
    ndata::Int
    function Data(x::T, y::T, err::T) where {T<:Vector{Real}}
        ndata = length(x)
        new(x, y, err, ndata)
    end
end

Data(df::DataFrame) = Data(df[:,1], df[:,2], df[:,3])
