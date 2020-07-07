__precompile__() # this module is safe to precompile
module IMinuit

export Minuit, migrad, minos, hesse
export func_argnames

include("Minuit.jl")

end
