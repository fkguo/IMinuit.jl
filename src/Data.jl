using DataFrames: DataFrame
using Plots
default(framestyle = :box, minorticks = true, fg_legend = cgrad(alpha = 0.3))

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
    plt_data(data::Data; xlab = "x", ylab = "y")

Make a errorbar plot of the data
"""
function plt_data1(data::Data; xlab = "x", ylab = "y", legend = :topleft)
    scatter(data.x, data.y, yerror = data.err, label = "Data", xlab = xlab, ylab = ylab, legend = legend)
end



"""
    chisq(dist::Function, data::Data, par; fitrange = ())

defines the χ² function: `fun` the function to be fitted to the data given by `data`.
The parameters are collected into `par`, given as an array or a tuple.

`fitrange`: default to the whole data set; may given as, e.g., `2:10`,
which means only fitting to the 2nd to the 10th data points
"""
function chisq(dist::Function, data::Data, par::Vector; fitrange = ())
    fitrange = (isempty(fitrange) && 1:data.ndata)
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
    plt_best(dist::Function, fit::Fit, data::Data; xrange = (), xlab = "x", ylab = "y", npt = 100)`

for plotting the comparison of the result from fit with the data.

`xrange`: range of `x` for plotting the best fit; if not given then use the range of `data.x`

`npt`: number of points computed for the best-fit curve, default = 100.
"""
function plt_best(dist::Function, fit::Fit, data::Data; xrange = (), xlab = "x", ylab = "y", npt = 100)
    paras1 = convert(Array, fit.args)

    dis(x) = (length(func_argnames(dist)) > 2 ? dist(x, paras1...) : dist(x, paras1))

    fig, ax = plt.subplots(figsize=(6,4))
    ax.minorticks_on(); ax.tick_params(which="both",direction="in", right="on", top="on")
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    xrange = (isempty(xrange) && data.x)
    wv = LinRange(xrange[1], xrange[end], npt)
    ax.errorbar(data.x, data.y, data.err, c = "C0", fmt="o", label = "Data")
    ax.plot(wv, dis.(wv), c = "C3", "-", label = "Best fit", lw=1.25)
    ax.legend();

    ax.set_xlim(wv[1], wv[end]); #ax.set_ylim(0, )
    ax.grid(true, alpha = 0.3)
end
