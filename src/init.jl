# IMinuit initialization

const minuit = PyNULL()
const mMinuit = PyNULL()
const iminuit_version = "1.4.5"

# initialization -- anything that depends on Python has to go here,
# so that it occurs at runtime (while the rest can be precompiled).
function __init__()
    copy!(minuit, pyimport_conda("iminuit", "iminuit=" * iminuit_version, "conda-forge"))
    copy!(mMinuit, pyimport(:iminuit).:Minuit)
end
