# IMinuit initialization

const iminuit = PyNULL()
const mMinuit = PyNULL()
const iminuit_version = "iminuit=2.18.0"
const cost = PyNULL()

# initialization -- anything that depends on Python has to go here,
# so that it occurs at runtime (while the rest can be precompiled).
# if PyCall is configured to use the Julia-specific Python, then iminuit
# can be automatically installed by pyimport_conda.
function __init__()
    copy!(iminuit, pyimport_conda("iminuit", "iminuit"))
    _version = iminuit.__version__
    if (_version < "2.0")
        println("The current iminuit version is " * _version * " which is not supported by this package, please upgrade to at least 2.0.0")
        println("update to the latest version?([Y]/n)")
        buf = readline()
        while buf != "y" && buf != "Y" && buf != ""
        	if buf == "n" || buf == "N"
                println("Selected not to update, aborting")
                exit()
            else
                println("Input \"$buf\" not supported")
                println("update to the latest version?([Y]/n)")
            end
            buf = readline()
        end
        run(`conda install iminuit`) # this updates the sys. conda, not the julia one.
        println("The iminuit has been successfully updated, quiting the process, please start again to refresh the environment")
        exit()
        #     copy!(iminuit, pyimport_conda("iminuit", iminuit_version, "conda-forge"))
        # #     _version = iminuit.:version.__version__
    end
    # println("iminuit version " * _version * " has been imported as `iminuit`." )
    copy!(mMinuit, pyimport_conda("iminuit", "iminuit").:Minuit)
    copy!(cost, pyimport_conda("iminuit.cost", "iminuit"))
    # The following converts the ArgsView type to Vector{Float64}, but takes too much time
    # it takes more than 50 μs
    # pytype_mapping(minuit._libiminuit.ArgsView, Vector{Float64})
    println("initiated!!")
end
