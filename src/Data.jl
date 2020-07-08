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
    function Data(x::Vector{T}, y::Vector{T}, err::Vector{T}) where {T<:Real}
        ndata = length(x)
        new(x, y, err, ndata)
    end
end

Data(df::DataFrame) = Data(df[:,1], df[:,2], df[:,3])

"""
    plt_data(data::Data; xlab = "x", ylab = "y")

Make a errorbar plot of the data
"""
function plt_data(data::Data; xlab = "x", ylab = "y")
    subplots(figsize = (6, 4))
    errorbar(data.x, data.y, data.err, fmt = "o", label = "Data" )
    legend()
    xlabel(xlab); ylabel(ylab)
    minorticks_on(); tick_params(which="both",direction="in", right="on", top="on")
    grid(which="major", axis="both", alpha=0.25)
end



"""
    χsq(fun, data::Data, par)

defines the χ² function: `fun` the function to be fitted to the data given by `data`.
The parameters are collected into `par`, given as an array or a tuple.
"""
function χsq(fun, data::Data, par; fitrange = 1:data.ndata)
    res = 0.0
    @simd for i = fitrange
        @inbounds res += ( (data.y[i]- fun(data.x[i], par))/data.err[i] )^2
    end
    return res
end
