function preprocess(fcn; kwds...)
    fitarg = ["error", "fix", "limit"]
    removed = ["errordef", "throw_nan", "print_level", "use_array_call"]
    if haskey(kwds, :name)
        args = Symbol.(kwds[:name])
    else
        args = func_argnames(fcn)
    end
    arg_dict = Dict(x => i for (i, x) in enumerate(args))
    len = length(args)
    stored_kwds = Dict{String,Any}()
    new_kwds = Vector()
    fitarg_dict = Dict{String,Array{Union{Float64,Bool,Tuple,Nothing}}}()
    fitarg_dict["error"] = fill(0.0, len)
    fitarg_dict["limit"] = fill((-Inf64, Inf64), len)
    fitarg_dict["fix"] = falses(len)
    # println("kwds ",kwds, "kwds")
    # println("kwds, $kwds")
    for (k, v) in kwds
        if v === nothing || k == :pedantic
            continue
        end
        k_str = String(k)
        if Symbol(k_str) in args
            push!(new_kwds, (k, v))
            continue
        elseif k_str in removed
            stored_kwds[k_str] = v
            continue
        elseif k_str in fitarg
            # println("values ", v)
            fitarg_dict[k_str] = v
            continue
        end
        udscore = findfirst('_', k_str)
        if udscore === nothing
            push!(new_kwds, (k, v))
            continue
        end
        typ = k_str[1:udscore-1]
        if typ in fitarg
            para = Symbol(k_str[udscore+1:end])
            fitarg_dict[typ][arg_dict[para]] = isa(v, Vector) ? Tuple(v) : v
        end
    end
    return (new_kwds, stored_kwds, fitarg_dict)
end

function preprocess(fcn, fit::AbstractFit; kwds...)
    fitarg = ["error", "fix", "limit"]
    removed = ["errordef", "throw_nan", "print_level", "use_array_call"]
    if haskey(kwds, :name)
        args = Symbol.(kwds[:name])
    else
        args = Symbol.(fit.parameters)
    end
    arg_dict = Dict(x => i for (i, x) in enumerate(args))
    stored_kwds = Dict{String,Any}()
    ini_value = Dict{Symbol, Float64}(args[i] => fit.values[i] for i in eachindex(args))
    new_kwds = Vector()
    fitarg_dict = Dict{String,AbstractVector{Union{Float64,Bool,Tuple,Nothing}}}("error" => [fit.errors...], "limit" => [fit.limits...], "fix" => [fit.fixed...])
    for (k, v) in kwds
        if v === nothing || k == :pedantic
            continue
        end
        k_str = String(k)
        if k_str in removed
            stored_kwds[k_str] = v
            continue
        elseif k_str in fitarg
            fitarg_dict[k_str] = v
            continue
        elseif k in args
            ini_value[k] = v
        end
        udscore = findfirst('_', k_str)
        if udscore === nothing
            push!(new_kwds, (k, v))
            continue
        end
        typ = k_str[1:udscore-1]
        if typ in fitarg
            para = Symbol(k_str[udscore+1:end])
            fitarg_dict[typ][arg_dict[para]] = isa(v, Vector) ? Tuple(v) : v
        end
    end
    if isa(fit, ArrayFit)
        ini_value = [ini_value[k] for k in args]
    end
    return (new_kwds, stored_kwds, ini_value, fitarg_dict)
end

