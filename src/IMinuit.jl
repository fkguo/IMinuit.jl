__precompile__() # this module is safe to precompile
module IMinuit

using PyCall
minuit = pyimport(:iminuit)

export Minuit, migrad, minos, hesse

include("Minuit.jl")

end
