using Combinatorics: combinations
using StatsBase: sample

@doc raw"""
    get_contours(fit::AbstractFit, χsq, parameters_combination::Vector{Int}; 
                 npts::Int=20, limits=true, sigma = 1.0)

For a given fit `fit` and the ``\chi^2`` function `χsq`, gives an array of parameter arrays,
with each array corresponding to a set of parameters obtained from calculating the `MINOS` ``1σ`` contour
(try to find `npts` points in the contour) for the two parameters in `parameters_combination`.

`parameters_combination` is an `Int` array of the numbers of that two parameters, e.g. it is `[1, 2]` for the first
two parameters and `[2, 3]` or the second and third parameters.

If `limits` is `true`, then fix one parameter to its bounds from `MINOS` of the best fit and get the values for
the other parameters; this runs over all parameters.
"""

# yes, you need to copy the state instead of set them in the init function
function copy_state!(from::T, to::T, cool_bits::Int) where {T<:AbstractFit}
    bitstream = bitstring(cool_bits)
    if bitstream[end] == '1'
        to.values = from.values
    end
    if bitstream[end-1] == '1'
        to.errors = from.errors
    end
    if bitstream[end-2] == '1'
        to.fixed = from.fixed
    end
    if bitstream[end-3] == '1'
        to.limits = from.limits
    end
end

function get_contours(fit::Fit, χsq, parameters_combination::Vector{Int}; npts::Int=20, limits=true, sigma=1.0)
    fit.valid || migrad(fit)        # if migrad has not been run then run it first
    length(fit.merrors) == 0 && minos(fit) # if minos has not been run then run it first
    fmin_1σ = fit.fval + 1.0   # χsq from the previous best fit + 1.0
    tol = 0.05  # tolerance allowing χ² to be ≤  fmin_1σ + tol
    # setting the parameters to the the best ones from the previous fit
    _fixed = fit.fixed
    narg = length(_fixed)
    _limits = [fit.limits...]
    _args = Symbol.(fit.parameters) #  func_argnames(χsq)
    @views for i in 1:length(_args)    # reset the limits of parameters to be the bounds from minos
        if !_fixed[i]
            _limits[i] = fit.values[i] .+ (fit.merrors[String(_args[i])][:lower], fit.merrors[String(_args[i])][:upper])
        end
    end

    container = Vector{Vector{Float64}}()
    ini_value = Dict{Symbol, Float64}(Symbol(_args[i]) => fit.values[i] for i in 1:narg)
    if limits
        @views for i in 1:narg
            if !_fixed[i]
                _fit1 = Minuit(χsq; ini_value...)
                _fit2 = Minuit(χsq; ini_value...)

                _fit1.values[i] = _limits[i][1]
                _fit2.values[i] = _limits[i][2]

                copy_state!(fit, _fit1, 6)
                copy_state!(fit, _fit2, 6)

                _fit1.fixed[i] = true
                _fit2.fixed[i] = true

                _fit1.limits = _limits
                _fit2.limits = _limits


                _fit1.strategy = 1
                migrad(_fit1)
                _fit2.strategy = 1
                migrad(_fit2)
                _fit1.fval < fmin_1σ + tol && push!(container, vcat(_fit1.fval, args(_fit1)))
                _fit2.fval < fmin_1σ + tol && push!(container, vcat(_fit2.fval, args(_fit2)))
            end
        end
    end

    # choose a pair of parameters, and try to compute the MINOS contour of npts points
    para1, para2 = _args[parameters_combination[1]], _args[parameters_combination[2]]
    contour_parameters = fit.mncontour(para1, para2, cl=sigma, size=npts)[1:end - 1, :]

    # for each contour point, get the values of the other parameters
    for pa in eachrow(contour_parameters)
        _fit1 = Minuit(χsq; ini_value...)
        _fit1.values[parameters_combination[1]], _fit1.values[parameters_combination[2]] = pa
        copy_state!(fit, _fit1, 6)
        _fit1.fixed[parameters_combination[1]] = true
        _fit1.fixed[parameters_combination[2]] = true


        #copy state
        _fit1.limits = _limits


        _fit1.strategy = 1
        migrad(_fit1)
        # filtering only those with χ² in 1σ
        _fit1.fval < fmin_1σ + tol && push!(container, vcat(_fit1.fval, args(_fit1)))
    end

    return container
end


function get_contours(fit::ArrayFit, χsq, parameters_combination::Vector{Int}; npts::Int=20, limits=true, sigma=1.0)
    fit.valid || migrad(fit)        # if migrad has not been run then run it first
    if length(fit.merrors) == 0 # if minos has not been run then run it first
        minos(fit)
    end
    fmin_1σ = fit.fval + 1.0   # χsq from the previous best fit + 1.0
    tol = 0.05  # tolerance allowing χ² to be ≤  fmin_1σ + tol
    # setting the parameters to the the best ones from the previous fit
    start_array = args(fit) #PyVector{Real}(fit.args) # does not support assignment
    _fixed = fit.fixed
    nargs = length(_fixed)
    _limits = [fit.limits...]
    _args = Symbol.(fit.parameters) # get names of parameters
    @views for i in 1:nargs    # reset the limits of parameters to be the bounds from minos
        if !_fixed[i]
            _limits[i] = fit.values[i] .+ (fit.merrors[String(_args[i])][:lower], fit.merrors[String(_args[i])][:upper])
        end
    end

    # igrad ? gradfun(par) = gradient(χsq, par) : gradfun = nothing

    container = Vector{Vector{Float64}}()
    # consider the upper and lower bounds for each parameter
    if limits
        @views for i in 1:nargs # for a in _args
            if !_fixed[i]
                start_array[i] = _limits[i][1]
                _fit1 = Minuit(χsq, start_array; name=_args)
                start_array[i] = _limits[i][2]
                _fit2 = Minuit(χsq, start_array; name=_args)

                copy_state!(fit, _fit1, 6)
                copy_state!(fit, _fit2, 6)

                _fit1.fixed[i] = true
                _fit2.fixed[i] = true

                #copy state
                _fit1.limits = _limits
                _fit2.limits = _limits

                _fit1.strategy = 1
                migrad(_fit1)
                _fit2.strategy = 1
                migrad(_fit2)
                _fit1.fval < fmin_1σ + tol && push!(container, vcat(_fit1.fval, args(_fit1)))
                _fit2.fval < fmin_1σ + tol && push!(container, vcat(_fit2.fval, args(_fit2)))
            end
        end
    end

    # choose a pair of parameters, and try to compute the MINOS contour of npts points
    para1, para2 = _args[parameters_combination[1]], _args[parameters_combination[2]]
    # println("*** ",fit, " ***")
    contour_parameters = fit.mncontour(para1, para2, cl=sigma, size=npts)[1:end - 1, :]

    # for each contour point, get the values of the other parameters
    for pa in eachrow(contour_parameters)
        start_array[parameters_combination[1]], start_array[parameters_combination[2]] = pa
        _fit1 = Minuit(χsq, start_array; name=_args)
        copy_state!(fit, _fit1, 6)
        _fit1.fixed[parameters_combination[1]] = true
        _fit1.fixed[parameters_combination[2]] = true

        _fit1.limits = _limits

        _fit1.strategy = 1
        migrad(_fit1)
        # println("*** ", _fit1.fval, " ***")
        # filtering only those with χ² in 1σ
        _fit1.fval < fmin_1σ + tol && push!(container, vcat(_fit1.fval, args(_fit1)))
    end

    return container
end


@doc raw"""
    get_contours_all(fit::AbstractFit, χsq; npts=20, limits=true, sigma = 1.0)

For a given fit `fit` and the ``\chi^2`` function `χsq`, gives a list of parameters sets which are at the edge of
``1σ`` `MINOS` contours for all combinations of varying parameters. The case of `limits` being `true` runs only once.
"""
function get_contours_all(fit::AbstractFit, χsq; npts=20, limits=true, sigma=1.0)
    fit.valid || migrad(fit)
    container = Vector()
    push!(container, vcat(fit.fval, args(fit))) # the first row set to tbe best fit
    free_pars_indices = Vector()
    for i in 1:length(fit.fixed)
        if fit.fixed[i] == false
            push!(free_pars_indices, i)
        end
    end
    for arr in combinations(free_pars_indices, 2)  # (1:npara, 2)
        push!(container, get_contours(fit, χsq, arr, npts=npts, limits=limits, sigma=sigma)...)
        limits = false  # run limits only once
    end
    # error("testing")
    return container
end


"""
    contour_df(fit::AbstractFit, χsq; npts=20, limits=true, sigma = 1.0)

parameters in the form of a dataframe.
"""
function contour_df(fit::AbstractFit, χsq; npts=20, limits=true, sigma=1.0)
    parameters_1sigma = get_contours_all(fit, χsq, npts=npts, limits=limits, sigma=sigma)
    # println(parameters_1sigma)
    df0 = DataFrame(parameters_1sigma, :auto)
    argnames = Array{Symbol,1}(undef, fit.npar)
    @. argnames = Symbol(fit.parameters)
    try
        DataFrame(collect.(eachrow(df0)), vcat(:chisq, argnames))
    catch
        @warn "No parameter sets were found."
    end
end


"""
    get_contours_given_parameter(fit::AbstractFit, χsq, para::T, range) 
        where {T <: Union{Symbol, String}}

gives parameter sets in one sigma for a given parameter constrained in a range.
If `fit` is an `ArrayFit` and no user-defined names have been given to the parameters, 
then `para` is `"x0"` or `:x0` for the 1st parameter, `"x1"` or `:x1` for the 2nd parameter, ...
"""
function get_contours_given_parameter(fit::Fit, χsq, para::T, range) where {T<:Union{Symbol,String}}
    fit.valid || migrad(fit)        # if migrad has not been run then run it first
    length(fit.merrors) == 0 && minos(fit) # if minos has not been run then run it first
    fmin_1σ = fit.fval + 1.0   # χsq from the previous best fit + 1.0
    # setting the parameters to the the best ones from the previous fit
    _args = Symbol.(fit.parameters) #  func_argnames(χsq)
    _fixed = fit.fixed
    nargs = length(_args)
    _limits = [fit.limits...]
    @views for i in 1:length(_args)    # reset the limits of parameters to be the bounds from minos
        if !_fixed[i]
            _limits[i] = fit.values[i] .+
                         (fit.merrors[String(_args[i])][:lower], fit.merrors[String(_args[i])][:upper])
        end
    end

    _para = findfirst(x -> x == String(para), fit.parameters)    # only the index would work in the new version, quite annoying
    container = Vector{Vector{Float64}}()
    ini_value = (Symbol(_args[i]) => fit.values[i] for i in 1::length(_args))
    for a in range
        _fit1 = Minuit(χsq; ini_value...)
        copy_state!(fit, _fit1, 6)
        _fit1.fixed[_para] = true
        _fit1.values[_para] = a

        _fit1.limits = _limits


        _fit1.strategy = 2
        _fit1.migrad()
        _fit1.fval ≤ fmin_1σ && push!(container, vcat(_fit1.fval, args(_fit1)))
    end

    return container
end


function get_contours_given_parameter(fit::ArrayFit, χsq, para::T, range) where {T<:Union{Symbol,String}}
    fit.valid || migrad(fit)        # if migrad has not been run then run it first
    length(fit.merrors) == 0 && minos(fit) # if minos has not been run then run it first
    fmin_1σ = fit.fval + 1.0   # χsq from the previous best fit + 1.0
    # setting the parameters to the the best ones from the previous fit
    start_array = args(fit) #PyVector{Real}(fit.args) # does not support assignment
    _args = Symbol.(fit.parameters)
    _fixed = fit.fixed
    nargs = length(_args)
    _limits = [fit.limits...]
    @views for i in 1:nargs    # reset the limits of parameters to be the bounds from minos
        if !_fixed[i]
            _limits[i] = fit.values[i] .+
                         (fit.merrors[String(_args[i])][:lower], fit.merrors[String(_args[i])][:upper])
        end
    end

    container = Vector{Vector{Float64}}()
    para_idx = findfirst(x -> x == Symbol(para), _args)
    for a in range
        start_array[para_idx] = a
        _fit1 = Minuit(χsq, start_array; name=_args)

        copy_state!(fit, _fit1, 6)
        _fit1.fixed[para_idx] = true
        _fit1.limits = _limits


        _fit1.strategy = 2
        migrad(_fit1)
        _fit1.fval ≤ fmin_1σ && push!(container, vcat(_fit1.fval, args(_fit1)))
    end

    return container
end


@doc """
    contour_df_given_parameter(fit::AbstractFit, χsq, para::T, range; limits = true) 
        where {T <: Union{Symbol, String}}

return parameter sets in one sigma for a given parameter constrained in a range as a `DataFrame`.

See also [`get_contours_given_parameter`](@ref).
"""
function contour_df_given_parameter(fit::AbstractFit, χsq, para::T, range; limits=true) where {T<:Union{Symbol,String}}
    parameters = get_contours_given_parameter(fit, χsq, para, range)
    df0 = DataFrame(parameters)
    argnames = Array{Symbol,1}(undef, fit.narg)
    @. argnames = Symbol(fit.parameters)
    try
        DataFrame(collect.(eachrow(df0)), vcat(:chisq, argnames))
    catch
        @warn "No parameter sets were found in the given range of $para. Try another range."
    end
end


@doc raw"""
    get_contours_samples(fit::AbstractFit, χsq, paras, ranges; nsamples = 100, MNbounds = true)

return 1σ parameter sets as an `Array` (the latter returns a `DataFrame`) for given parameters constrained in `ranges`:
* if `paras` is a single parameter, then take equally spaced `nsamples` in `ranges` given in the form of `(min, max)`;
* if `paras` contain more (``\geq 2``) parameters, then `paras` should be of the form `(:para1, :para2)`, `ranges` should be of the form `((min1, max1), (min2, max2))`,
and values for the parameters given in `paras` are randomly sampled in the given `ranges`;
* if `MNbounds` is true, then constrain the parameters in the range provided by `MINOS` no matter whether that is valid or not (to be improved by checking the validity)
* if `igrad` is true, then use `ForwardDiff.gradient` to compute the gradient.
* For using array parameters, if no user-defined names have been given to the parameters,
`paras` should be given such that `"x0"` or `:x0` for the 1st parameter, `"x1"` or `:x1` for the 2nd parameter, ...
"""
function get_contours_samples(fit::Fit, χsq, paras, ranges; nsamples=100, MNbounds=true, igrad=false)
    fit.valid || migrad(fit)        # if migrad has not been run then run it first
    length(fit.merrors) == 0 && minos(fit) # if minos has not been run then run it first
    fmin_1σ = fit.fval + 1.0   # χsq from the previous best fit + 1.0
    # setting the parameters to the the best ones from the previous fit


    _args = Symbol.(fit.parameters)
    _fixed = fit.fixed
    nargs = length(_args)
    _limits = [fit.limits...]
    ini_value = (_args[i] => fit.values[i] for i in 1:nargs)
    if MNbounds
        @views for i in 1:length(_args)    # reset the limits of parameters to be the bounds from minos
            if !_fixed[i]
                _limits[i] = fit.values[i] .+
                             (fit.merrors[String(_args[i])][:lower], fit.merrors[String(_args[i])][:upper])
            end
        end
    end

    container = Vector{Vector{Float64}}()
    push!(container, vcat(fit.fval, args(fit))) # the first row set to tbe best fit

    # make range for each parameter with length len; if one only 1 parameter
    @views if isa(ranges[1], Number)
        sam = LinRange(ranges[1], ranges[end], nsamples)
    else
        r = map(x -> LinRange(x[1], x[end], nsamples), ranges)
        sam = zip(sample.(r, nsamples)...)
    end

    igrad ? gradf(a...) = gradient(x -> χsq(x...), [a...]) : gradf = nothing

    par_to_idx = Dict(x => i for (i, x) in enumerate(fit.parameters))


    @inbounds @simd for p in unique(sam)
        @views if isa(ranges[1], Number)
            _para = findfirst(x -> x == paras, _args)    # only the index would work in the new version, quite annoying
            _fit1 = Minuit(χsq; grad=gradf, ini_value...)
            copy_state!(fit, _fit1, 6)
            _fit1.values[_para] = p
            _fit1.fixed[_para] = true
            _fit1.limits = _limits
        else
            _fit1 = Minuit(χsq; grad=gradf, ini_value...)
            for (idx, para) in enumerate(paras)
                _fit1.fixed[par_to_idx[String(para)]] = true
                _fit1.values[par_to_idx[String(para)]] = p[idx]
            end
        end
        _fit1.strategy = 1
        migrad(_fit1)
        _fit1.fval ≤ fmin_1σ && push!(container, vcat(_fit1.fval, args(_fit1)))
    end

    return container
end


function get_contours_samples(fit::ArrayFit, χsq, paras, ranges; nsamples=100, MNbounds=true, igrad=false)
    fit.valid || migrad(fit)        # if migrad has not been run then run it first
    length(fit.merrors) == 0 && minos(fit) # if minos has not been run then run it first
    fmin_1σ = fit.fval + 1.0   # χsq from the previous best fit + 1.0
    start_array = args(fit)
    _args = Symbol.(fit.parameters)
    if _args[1] != :x0
        if isa(paras, Tuple) || isa(paras, Vector)
            let para = String(paras[1])
            	if para[1] == 'x' && isascii(para[end])
                	paras = collect(_args[parse(Int, String(paras[i])[2:end]) + 1] for i in eachindex(paras))
              end 
            end
        elseif String(paras)[1] == 'x' && isascii(String(paras)[end])
            paras = _args[parse(Int, String(paras)[2:end]) + 1]
        end
    end
    # println(_args)
    _fixed = fit.fixed
    nargs = length(_args)
    _limits = [fit.limits...]
    if MNbounds
        @views for i in 1:length(_args)    # reset the limits of parameters to be the bounds from minos
            if !_fixed[i]
                _limits[i] = fit.values[i] .+
                             (fit.merrors[String(_args[i])][:lower], fit.merrors[String(_args[i])][:upper])
            end
        end
    end

    container = Vector{Vector{Float64}}()
    push!(container, vcat(fit.fval, args(fit))) # the first row set to tbe best fit

    # make range for each parameter with length len; if one only 1 parameter
    @views if isa(ranges[1], Number)
        sam = LinRange(ranges[1], ranges[2], nsamples)
        pfix = Symbol(:fix_, paras) => true
    else
        r = map(x -> LinRange(x[1], x[2], nsamples), ranges)
        sam = zip(sample.(r, nsamples)...)
    end

    igrad ? gradf(a...) = gradient(x -> χsq(x...), [a...]) : gradf = nothing

    # if using @inbounds, then crashes julia
    for p in unique(sam)
        @views if isa(ranges[1], Number)
            _para = findfirst(x -> x == paras, _args)    # only the index would work in the new version, quite annoying
            start_array[_para] = p
            _fit1 = Minuit(χsq, start_array; name=_args, grad=gradf)
            copy_state!(fit, _fit1, 6)
            _fit1.fixed[_para] = true
            _fit1.limits = _limits
        else
            _arg_to_idx = Dict(x => i for (i, x) in enumerate(_args))
            # par_to_idx = Dict(x => i for (i, x) in enumerate(fit.parameters))
            @views for i in eachindex(paras)
                start_array[_arg_to_idx[paras[i]]] = p[i]
            end
            _fit1 = Minuit(χsq, start_array; name=_args, grad=gradf)
            copy_state!(fit, _fit1, 6)
            _fit1.limits = _limits
            for para in paras
                _fit1.fixed[_arg_to_idx[para]] = true
            end
        end
        _fit1.strategy = 1
        migrad(_fit1)
        _fit1.fval ≤ fmin_1σ && push!(container, vcat(_fit1.fval, args(_fit1)))
    end

    return container
end

@doc raw"""
    contour_df_samples(fit::AbstractFit, χsq, paras, ranges; nsamples = 100, MNbounds=true)

gives  1σ parameter sets as a `DataFrame` for given parameters constrained in `ranges`:
* if `paras` is a single parameter, then take equally spaced `nsamples` in `ranges` given in the form of `(min, max)`;
* if `paras` contain more parameters, then `paras` should be of the form `(:para1, :para2)`, `ranges` should be of the form `((min1, max1), (min2, max2))`;
* `paras` can be more than 2. Values for the parameters given in `paras` are randomly sampled in the given `ranges`.
* is `MNbounds` is true, then constrain the parameters in the range provided by `MINOS` no matter whether that is valid or not (to be improved by checking the validity)
* if `igrad` is true, then use `ForwardDiff.gradient` to compute the gradient
"""
function contour_df_samples(fit::AbstractFit, χsq, paras, ranges; nsamples=100, MNbounds=true, igrad=false)
    # println("paras! $paras")
    parameters = get_contours_samples(fit, χsq, paras, ranges; nsamples=nsamples, MNbounds=MNbounds, igrad=igrad)
    df0 = DataFrame(parameters, :auto)
    argnames = Array{Symbol,1}(undef, fit.npar)
    @. argnames = Symbol(fit.parameters)
    try
        DataFrame(collect.(eachrow(df0)), vcat(:chisq, argnames))
    catch
        @warn "No parameter sets were found in the given ranges of $paras. Try another ranges."
    end
end
