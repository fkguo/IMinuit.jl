# IMinuit initialization

const minuit = PyNULL()
const mMinuit = PyNULL()

# initialization -- anything that depends on Python has to go here,
# so that it occurs at runtime (while the rest can be precompiled).
function __init__()
    copy!(minuit, pyimport_conda("iminuit", "iminuit", "conda-forge"))
    copy!(mMinuit, pyimport(:iminuit).:Minuit)
end
