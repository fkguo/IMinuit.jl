# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.7
#   kernelspec:
#     display_name: Julia 1.9.1
#     language: julia
#     name: julia-1.9
# ---

# +

using Pkg
Pkg.develop(path = "/home/akizuki/IMinuit.jl")

# -

using IMinuit
using BenchmarkTools

f(x) = x[1]^2 + (x[2]-1)^2 + (x[3]-2)^4
f1(x, y, z) = x^2 + (y-1)^2 + (z-2)^4

# using array parameters
m = Minuit(f, [1, 1, 4], fix_x0 = true)
migrad(m)

# a new fit can continue from the previous fit
m_new = Minuit(f, m, fix_x0 = false)
migrad(m_new)

# using array parameters, using `ForwardDiff: gradient` to compute the gradient
gradf(x) = gradient(f, x)
mgrad = Minuit(f, [1, 1, 4], grad = gradf)
migrad(mgrad)

# parameters are given individually
m1 = Minuit(f1, x = 1, y = 1, z = 4)
migrad(m1)

@show iminuit.__version__
propertynames(m)

# the doc strings are from `iminuit`
@doc migrad

# ## Example: Fit to the BES data of the π⁺π⁻ energy distribution of ψ'→J/ψπ⁺π⁻
#
# The data are taken from [BES Collaboration, Phys. Rev. D 62 (2000) 032002](https://inspirehep.net/literature/507637).
#
# Here we use a simple model, which is not meant to be the correct one, to fit to the data.



using CSV
using DataFrames
using Plots
pyplot(framestyle = :box, minorticks = 5)
using LaTeXStrings

data_df = DataFrame(CSV.File("./testdata.csv"))
const data = Data(data_df)

@plt_data data  xlab=L"m_{\pi\pi}"*" [GeV]"  ylab="Events"

# +
const M = 3.686; const mπ = 0.14; const mJ = 3.097; 

λ(x, y, z) = x^2 + y^2 + z^2 - 2x*y - 2y*z - 2z*x

# a simple function that will be used to fit the data: QCD multipole expansion model for ψ'→J/ψπ⁺π⁻
# The important ππ FSI effect is not taken into account
# bg is just for introducing a third parameter
function dist(w, N, c, bg) 
    if (w ≤ 2mπ || w ≥ M-mJ)
        res = 0.0
    else
        q1 = sqrt(λ(w^2, mπ^2, mπ^2))/(2w)
        q2 = sqrt(λ(M^2, w^2, mJ^2))/(2M)
        res = N * q1 * q2 * (w^2 - c*mπ^2)^2 + bg
    end
    return res * 1e6
end;

dist(x, p) = dist(x, p...)

# +
# parameters given individually
χsq1(N, c, bg) = chisq(dist, data, (N, c, bg));

# all parameters are vairables of χsq
fit1 = Minuit(χsq1, N = 1, c = 2, bg = 0, error_N = 0.1, error_c = 0.1, error_bg = 0.1)
fit1.strategy = 1;

# +
# parameters are collected into a tuple or an array, which is the only variable of χsq
parname = [:N, :c, :bg]
χsq(par) = chisq(dist, data, par)
gradf(par) = gradient(χsq, par)
fit = Minuit(χsq, [1, 2, 0], error = 0.1*ones(3), name = parname, grad = gradf)
fit.strategy = 1;

# or simply using model_fit or @model_fit
fit2 = model_fit(dist, data, [1, 2, 0], error = 0.1*ones(3), name = parname, fix_bg = true)
fit2.strategy = 1;
migrad(fit2)
# -

# the privous fit status can be passed to a new fit
fit2_new = model_fit(dist, data, fit2, name = parname, fix_bg = false)
migrad(fit2_new)

@btime migrad(fit)
minos(fit)
migrad(fit)

@btime migrad(fit1)
minos(fit1)
migrad(fit1)

minos(fit1)

# the ordering of dist, fit and data does not matter
@plt_best dist fit data

# MIGRAD contour of two parameters with the other ones fixed
# needs PyPlot or the pyplot backend of Plots
draw_contour(fit1,:N, :c, bound=3, bins=100)

# contour of parameter space from MINOS
draw_mncontour(fit1,:N, :c, nsigma=3, numpoints=100);draw_mncontour(fit1,:N, :c, nsigma=2, numpoints=100)
draw_mncontour(fit1,:N, :c, nsigma=1, numpoints=100)

matrix(fit1, correlation = true)

@show fit1.matrix
matrix(fit1)

# this gives parameter sets at the 1σ boundary
@time contour_df(fit1, χsq1, npts = 5)

# random sampling of parameters in given ranges, keeping those within 1σ
@time parsam_df = contour_df_samples(fit1, χsq1, (:N, :c), ([2.5,2.8], [4.0,4.3]), nsamples = 3000)

# get parameter ranges
extrema(parsam_df.:N), extrema(parsam_df.:c)

scatter(parsam_df.:N, parsam_df.:c, xlab = "N", ylab = "c")

@time contour_df_samples(fit, χsq, (:x0, :x1, :x2), ([2.5,2.8],  [4.0,4.3], (-3e-5,3e-5)), nsamples = 1000)

@time contour_df_samples(fit, χsq, :x0, (2.5,2.8), nsamples = 20)

# ## Example A in the iminuit tutorial
#
# The example is [Example A: Fit of a gaussian model to a histogram](https://nbviewer.jupyter.org/github/scikit-hep/iminuit/blob/master/tutorial/automatic_differentiation.ipynb)



using PyCall

# import numpy from Python to generate the same data as in the example
np = pyimport(:numpy)
default_rng = pyimport("numpy.random").default_rng
rng = default_rng(seed=1)
const w, xe = np.histogram(rng.normal(0, 1, 10000), bins=1000)

# +
# define the model and the score function to minimize
using SpecialFunctions

function cdf(x, par)
    mu, sigma = par
    z = (x - mu) / sigma
    return 0.5 * (1 + erf(z / sqrt(2))) 
end

function score(par)
    amp = par[1]
    rest = par[2:end]
    mu = amp * (cdf.(xe[2:end], Ref(rest)) - cdf.(xe[1:end-1], Ref(rest)) )
    return 2 * sum(@. mu - w * log(mu + 1e-100))
end
# -

const start_values = [1.5 * sum(w), 1.0, 2.0]
const limits = [(0, nothing), nothing, (0, nothing)];

# +
# w/o grad
m = Minuit(score, start_values, limit=limits)
m.strategy = 0

# using grad
grad_fd(pars) =  gradient(score, pars)
m_fd = Minuit(score, start_values, limit=limits, grad = grad_fd)
m_fd.strategy = 0;
# -

@btime migrad(m)

@btime migrad(m_fd)


