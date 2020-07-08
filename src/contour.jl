using Combinatorics
import StatsBase: sample

"""
    get_contours(fit::Fit, χsq, parameters_combination::Vector{Int}; npts::Int=20, limits=true, sigma = 1.0)

For a given fit `fit` and the \$\\chi^2\$ function `χsq`, gives an array of parameter arrays,
with each array corresponding to a set of parameters obtained from calculating the `MINOS` \$1\\sigma\$ contour
(try to find `npts` points in the contour) for the two parameters in `parameters_combination`.

`parameters_combination` is an `Int` array of the numbers of that two parameters, e.g. it is `[1, 2]` for the first
two parameters and `[2, 3]` or the second and third parameters.

If `limits` is `true`, then fix one parameter to its bounds from `MINOS` of the best fit and get the values for
the other parameters; this runs over all parameters.
"""
function get_contours(fit::Fit, χsq, parameters_combination::Vector{Int}; npts::Int=20, limits=true, sigma = 1.0)
    fit.migrad_ok() || fit.migrad();        # if migrad has not been run then run it first
    length(fit.merrors) == 0 && fit.minos() # if minos has not been run then run it first
    fmin_1σ = fit.fval + 1.0;   # χsq from the previous best fit + 1.0
    tol = 0.05  # tolerance allowing χ² to be ≤  fmin_1σ + tol
    # setting the parameters to the the best ones from the previous fit
    kwdarg = Dict{Symbol, Any}(Symbol(k) => v for (k,v) in fit.fitarg)
    args = func_argnames(χsq)
    @views for i in 1: length(args)    # reset the limits of parameters to be the bounds from minos
        kwdarg[Symbol(:limit_, args[i])] = get(fit.args, i-1) .+
             (fit.merrors[String(args[i]), -1], fit.merrors[String(args[i]), 1])
    end

    container = Vector{Vector{Float64}}()
    # consider the upper and lower bounds for each parameter
    if limits
        @views for a in args
            fit1 = Minuit(χsq; kwdarg..., a => kwdarg[Symbol(:limit_, a)][1], Symbol("fix_", a) => true)
            fit2 = Minuit(χsq; kwdarg..., a => kwdarg[Symbol(:limit_, a)][2], Symbol("fix_", a) => true)
            fit1.strategy = 2; fit1.migrad()
            fit2.strategy = 2; fit2.migrad()
            fit1.fval  < fmin_1σ + tol && push!(container, vcat(fit1.fval, convert(Array, fit1.values)) )
            fit2.fval  < fmin_1σ + tol && push!(container, vcat(fit2.fval, convert(Array, fit2.values)) )
        end
    end

    # choose a pair of parameters, and try to compute the MINOS contour of npts points
    @views para1, para2 = args[parameters_combination[1]], args[parameters_combination[2]]
    # somehow when np≥6, cannot get results for [2,4]...; but succeeded before, don't know why
#     npts = parameters_combination == [2,4] ? 5 : npts
    @views contour_parameters = fit.mncontour(para1, para2, sigma = sigma, numpoints=npts )[3]

    # for each contour point, get the values of the other parameters
    for pa in contour_parameters
        dict = Dict(zip((para1, para2), pa) ) # using dict to overwrite those in kwdarg
        fit1 = Minuit(χsq; kwdarg..., dict...,  #errordef=1,
                  Symbol("fix_", para1) => true, Symbol("fix_", para2) => true)
        fit1.strategy = 2; fit1.migrad()
        # filtering only those with χ² in 1σ
        fit1.fval  < fmin_1σ + tol && push!(container, vcat(fit1.fval, convert(Array, fit1.values)) )
    end

    return container
end


"""
    get_contours_all(fit::Fit, χsq; npts=20, limits=true, sigma = 1.0)

For a given fit `fit` and the \$\\chi^2\$ function `χsq`, gives a list of parameters sets which are at the edge of
\$1\\sigma\$ `MINOS` contours for all parameters combinations. The case of `limits` being `true` runs only once.
"""
function get_contours_all(fit::Fit, χsq; npts=20, limits=true, sigma = 1.0)
    npara = length(func_argnames(χsq))
    container = Vector()
    push!(container, vcat(fit.fval, convert(Array, fit.values)) ) # the first row set to tbe best fit
    for arr in combinations(1:npara, 2)
        push!(container, get_contours(fit, χsq, arr, npts=npts, limits=limits, sigma = sigma)...)
        limits = false  # run limits only once
    end
    return container
end


"""
    contour_df(fit::Fit, χsq; npts=20, limits=true, sigma = 1.0)

parameters in the form of a dataframe
"""
function contour_df(fit::Fit, χsq; npts=20, limits=true, sigma = 1.0)
    parameters_1sigma =  get_contours_all(fit, χsq, npts=npts, limits=limits, sigma = sigma)
    df0 = DataFrame(parameters_1sigma)
    return DataFrame( collect.(eachrow(df0)), vcat(:χ², func_argnames(χsq)) )
end


"""
    get_contours_given_parameter(fit::Fit, χsq, para::T, range) where {T <: Union{Symbol, String}}

gives parameter sets in one sigma for a given parameter constrained in a range
"""
function get_contours_given_parameter(fit::Fit, χsq, para::T, range) where {T <: Union{Symbol, String}}
    fit.migrad_ok() || fit.migrad();        # if migrad has not been run then run it first
    length(fit.merrors) == 0 && fit.minos() # if minos has not been run then run it first
    fmin_1σ = fit.fval + 1.0;   # χsq from the previous best fit + 1.0
    # setting the parameters to the the best ones from the previous fit
    kwdarg = Dict{Symbol, Any}(Symbol(k) => v for (k,v) in fit.fitarg)
    args = func_argnames(χsq)
    @views for i in 1: length(args)    # reset the limits of parameters to be the bounds from minos
        kwdarg[Symbol(:limit_, args[i])] = get(fit.args, i-1) .+
             (fit.merrors[String(args[i]), -1], fit.merrors[String(args[i]), 1])
    end

    container = Vector{Vector{Float64}}()
    for a in range
        fit1 = Minuit(χsq; kwdarg..., para => a, Symbol(:fix_, para) => true)
        fit1.strategy = 2; fit1.migrad()
        fit1.fval  ≤ fmin_1σ && push!(container, vcat(fit1.fval, convert(Array, fit1.values)) )
    end

    return container
end


"""
    contour_df_given_parameter(fit, χsq, para::T, range; limits = true) where {T <: Union{Symbol, String}}

gives parameter sets in one sigma for a given parameter constrained in a range as a `DataFrame`
"""
function contour_df_given_parameter(fit::Fit, χsq, para::T, range; limits = true) where {T <: Union{Symbol, String}}
    parameters =  get_contours_given_parameter(fit, χsq, para, range)
    df0 = DataFrame(parameters)
    return DataFrame( collect.(eachrow(df0)), vcat(:χ², func_argnames(χsq)) )
end


"""
    get_contours_samples(fit, χsq, paras, ranges; nsamples = 100, MNbounds = true)

gives \$1\\sigma\$ parameter sets as an `Array` for given parameters constrained in `ranges`:
* if `paras` is a single parameter, then take equally spaced `nsamples` in `ranges` given in the form of `(min, max)`;
* if `paras` contain more parameters, then `paras` should be of the form `(:para1, :para2)`, `ranges` should be of the form `((min1, max1), (min2, max2))`;
* `paras` can be more than 2. Values for the parameters given in `paras` are randomly sampled in the given `ranges`.
* if `MNbounds` is true, then constrain the parameters in the range provided by `MINOS` no matter whether that is valid or not (to be improved by checking the validity)
* if `igrad` is true, then use ForwardDiff.gradient to compute the gradient
"""
function get_contours_samples(fit::Fit, χsq, paras, ranges; nsamples = 100, MNbounds = true, igrad = false)
    fit.migrad_ok() || fit.migrad();        # if migrad has not been run then run it first
    length(fit.merrors) == 0 && fit.minos() # if minos has not been run then run it first
    fmin_1σ = fit.fval + 1.0;   # χsq from the previous best fit + 1.0
    # setting the parameters to the the best ones from the previous fit
    kwdarg = Dict{Symbol, Any}(Symbol(k) => v for (k,v) in fit.fitarg)
    args = func_argnames(χsq)
    if MNbounds
        @views for i in 1: length(args)    # reset the limits of parameters to be the bounds from minos
            kwdarg[Symbol(:limit_, args[i])] = get(fit.args, i-1) .+
                 (fit.merrors[String(args[i]), -1], fit.merrors[String(args[i]), 1])
        end
    end

    container = Vector{Vector{Float64}}()
    push!(container, vcat(fit.fval, convert(Array, fit.values)) ) # the first row set to tbe best fit

    # make range for each parameter with length len; if one only 1 parameter
    @views if isa(ranges[1], Number)
        sam = LinRange(ranges[1], ranges[2], nsamples)
        pfix = Symbol(:fix_, paras) => true
    else
        r = map(x -> LinRange(x[1], x[2], nsamples), ranges)
        sam = zip(sample.(r, nsamples)...)
        pfix = map( x -> Symbol(:fix_, x) => true, paras)
    end

    igrad ? gradf(a...) = gradient(x->χsq(x...), [a...]) : gradf = nothing

    @inbounds @simd for p in unique(sam)
        if isa(ranges[1], Number)
            fit1 = Minuit(χsq; kwdarg..., paras => p, Symbol(:fix_, paras) => true, grad = gradf)
        else
            pdict = Dict(zip(paras, p))
            fit1 = Minuit(χsq; kwdarg..., pdict..., pfix..., grad = gradf)
        end
        fit1.strategy = 2; fit1.migrad()
        fit1.fval  ≤ fmin_1σ && push!(container, vcat(fit1.fval, convert(Array, fit1.values)) )
    end

    return container
end

"""
    contour_df_samples(fit::Fit, χsq, paras, ranges; nsamples = 100, MNbounds=true)

gives \$1\\sigma\$ parameter sets as a `DataFrame` for given parameters constrained in `ranges`:
* if `paras` is a single parameter, then take equally spaced `nsamples` in `ranges` given in the form of `(min, max)`;
* if `paras` contain more parameters, then `paras` should be of the form `(:para1, :para2)`, `ranges` should be of the form `((min1, max1), (min2, max2))`;
* `paras` can be more than 2. Values for the parameters given in `paras` are randomly sampled in the given `ranges`.
* is `MNbounds` is true, then constrain the parameters in the range provided by `MINOS` no matter whether that is valid or not (to be improved by checking the validity)
* if `igrad` is true, then use ForwardDiff.gradient to compute the gradient
"""
function contour_df_samples(fit::Fit, χsq, paras, ranges; nsamples = 100, MNbounds=true, igrad = false)
    parameters = get_contours_samples(fit, χsq, paras, ranges, nsamples = nsamples, MNbounds = MNbounds, igrad = igrad)
    df0 = DataFrame(parameters)
    return DataFrame( collect.(eachrow(df0)), vcat(:χ², func_argnames(χsq)) )
end
