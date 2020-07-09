# IMinuit initialization

const minuit = PyNULL()
const mMinuit = PyNULL()
const iminuit_version = "1.4.5" # this has to be at least 1.4.1

# initialization -- anything that depends on Python has to go here,
# so that it occurs at runtime (while the rest can be precompiled).
# if PyCall is configured to use the Julia-specific Python, then iminuit
# can be automatically installed by pyimport_conda.
function __init__()
    copy!(minuit, pyimport_conda("iminuit", "iminuit=" * iminuit_version, "conda-forge"))
    copy!(mMinuit, pyimport(:iminuit).:Minuit)
end
