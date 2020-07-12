using DataFrames: DataFrame
# using Plots: scatter, plot!, default
# default(framestyle = :box, minorticks = 4)

"""
    Data(x::T, y::T, err::T) where {T<:Vector{Real}}
    Data(df::DataFrame)

Fields: `x, y, err, ndata`

This defines a type for data with three columns:` x, y, err`; `ndata` is the number of data rows.

Only symmetric errors (of `y`) are supported.
"""
struct Data
    x::Vector{Float64}
    y::Vector{Float64}
    err::Vector{Float64}
    ndata::Int
    function Data(x::Vector{T}, y::Vector{T}, err::Vector{T}) where {T<:Real}
        ndata = length(x)
        new(x, y, err, ndata)
    end
end

Data(df::DataFrame) = Data(df[:,1], df[:,2], df[:,3])

"""
    @plt_data(data, kws...)

A convenient mascro to make an errorbar plot of the `data`; all combinations of
keyword settings for `scatter` in `Plots` can be used for the optional arguments `kws...`
"""
macro plt_data(data, kws...)
    _plt = quote
        if isempty($kws)
            scatter($data.x, $data.y, yerror = $data.err,
            xlab = "x", ylab = "y", label = "Data")
        else
            scatter($data.x, $data.y, yerror = $data.err; $(kws...) )
        end
    end
    return esc(_plt)
end



"""
    chisq(dist::Function, data::Data, par; fitrange = ())

defines the χ² function: `fun` the function to be fitted to the data given by `data`.
The parameters are collected into `par`, given as an array or a tuple.

`fitrange`: default to the whole data set; may given as, e.g., `2:10`,
which means only fitting to the 2nd to the 10th data points.
"""
function chisq(dist::Function, data::Data, par::AbstractVector; fitrange = ())
    fitrange = (isempty(fitrange) && 1:data.ndata)
    # __dist__(x, par) = (length(func_argnames(dist)) > 2 ? dist(x, par...) : dist(x, par) )
    res = 0.0
    @simd for i = fitrange
        @inbounds res += ( (data.y[i]- dist(data.x[i], par))/data.err[i] )^2
    end
    return res
end
function chisq(dist::Function, data::Data, par::Tuple; fitrange = ())
    fitrange = (isempty(fitrange) && 1:data.ndata)
    res = 0.0
    @simd for i = fitrange
        @inbounds res += ( (data.y[i]- dist(data.x[i], par))/data.err[i] )^2
    end
    return res
end

"""
    @plt_best(dist, fit, data, kws...)

A convenient macro for comparing the best-fit result with the data; all combinations of
keyword settings for `plot` in `Plots` can be used for the optional arguments `kws...`
"""
macro plt_best(dist, fit, data, kws...)
    _expr = quote
        _npts = 100
        _paras = args($fit)
        _dis(x) = (typeof(fit) == ArrayFit ? $dist(x, _paras) : $dist(x, _paras...) )
        _xrange = $data.x #(isempty($xrange) ? $data.x : $xrange)
        _wv = LinRange(_xrange[1], _xrange[end], _npts)

        scatter($data.x, $data.y, yerror = $data.err, label = "Data")
        if isempty($kws)
            plot!(_wv, _dis.(_wv), xlab = "x", ylab = "y", label = "Best fit", lw = 1.5 )
        else
            plot!(_wv, _dis.(_wv); $(kws...)  )
        end
    end
    esc( _expr )
end
