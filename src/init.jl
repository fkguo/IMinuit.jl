# IMinuit initialization

const iminuit = PyNULL()
const mMinuit = PyNULL()
const iminuit_version = "iminuit" # "iminuit=1.4.9"

# initialization -- anything that depends on Python has to go here,
# so that it occurs at runtime (while the rest can be precompiled).
# if PyCall is configured to use the Julia-specific Python, then iminuit
# can be automatically installed by pyimport_conda.
function __init__()
    copy!(iminuit, pyimport_conda("iminuit", iminuit_version, "conda-forge"))
    # _version = iminuit.:version.__version__
    # if (_version < "1.4.1") | (_version > "1.4.9")
    #     # println("iminuit version" * _version * "1.4.1. It will be updated now.")
    #     println("iminuit v1.4.9 will be installed now.")
    #     run(`conda install iminuit=1.4.9 -c conda-forge`) # this updates the sys. conda, not the julia one.
    #     copy!(iminuit, pyimport_conda("iminuit", iminuit_version, "conda-forge"))
    # #     _version = iminuit.:version.__version__
    # end
    # println("iminuit version " * _version * " has been imported as `iminuit`." )
    copy!(mMinuit, pyimport_conda("iminuit", iminuit_version, "conda-forge").:Minuit)
    # The following converts the ArgsView type to Vector{Float64}, but takes too much time
    # it takes more than 50 μs
    # pytype_mapping(minuit._libiminuit.ArgsView, Vector{Float64})
end
