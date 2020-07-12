# IMinuit initialization

const minuit = PyNULL()
const mMinuit = PyNULL()
const iminuit_version = "iminuit=1.4.6" # this has to be at least 1.4.1

# initialization -- anything that depends on Python has to go here,
# so that it occurs at runtime (while the rest can be precompiled).
# if PyCall is configured to use the Julia-specific Python, then iminuit
# can be automatically installed by pyimport_conda.
function __init__()
    copy!(minuit, pyimport_conda("iminuit", iminuit_version, "conda-forge"))
    copy!(mMinuit, pyimport_conda("iminuit", iminuit_version, "conda-forge").:Minuit)
    # The following converts the ArgsView type to Vector{Float64}, but takes too much time
    # it takes more than 50 μs
    # pytype_mapping(minuit._libiminuit.ArgsView, Vector{Float64})
end
