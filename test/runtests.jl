using IMinuit
using Test
using DataFrames, CSV

# @testset "IMinuit.jl" begin
    # Write your tests here.
# end

@testset "minimizing simple functions: " begin
    f(x) = x[1]^2 + (x[2]-1)^2 + (x[3]-2)^4
    f1(x, y, z) = x^2 + (y-1)^2 + (z-2)^4

    m = Minuit(f, [1,1,0])
    m1 = Minuit(f1, x = 1, y = 1, z = 0)

    migrad(m)
    migrad(m1)

    @test typeof(m) == ArrayFit
    @test typeof(m1) == Fit

    # the second parameter should be about 1
    @test get(m.values, 1) ≈ 1
    @test get(m1.values, 1) ≈ 1

    @test length(get_contours_all(m, f, npts = 3)) == 19
    @test length(get_contours_all(m1, f1, npts = 3)) == 19
end

@testset "a simple but realistic fit: " begin
    datadf = DataFrame!(CSV.File("../docs/testdata.csv"))
    data = Data(datadf)

    M = 3.686; mπ = 0.14; mJ = 3.097;

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
    end

    # dist(w, par) = dist(w, par...);

    χsq(par) = chisq(dist, data, par)
    χsq1(N, c, bg) = chisq(dist, data, (N, c, bg));

    gradf(par) = gradient(χsq, par)
    fit = Minuit(χsq, [1, 2, 0], error = 0.1*ones(3), grad = gradf)
    fit.strategy = 0;

    fit1 = Minuit(χsq1, N = 1, c = 2, bg = 0, error_N = 0.1, error_c = 0.1, error_bg = 0.1)
    fit1.strategy = 0;

    migrad(fit)

    migrad(fit1)

    @test isapprox(get(fit.values, 0), 2.61, rtol = 1e-2)
    @test isapprox(get(fit1.values, 0), 2.61, rtol = 1e-2)

end
