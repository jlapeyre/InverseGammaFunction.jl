using SpecialFunctions: loggamma
using InverseGammaFunction
using Test


@testset "InverseGammaFunction.jl" begin
    @test isapprox(invgamma(100.0), 5.892518696343773)
    @test isapprox(invgamma(100.0; branch=:lo), 0.00994357323626315)

    for y in (1.0, 2.0, 10.0, 100.0, 0.0, 0.1)
        @test isapprox(loggamma(invloggamma(y)), y)
        @test isapprox(loggamma(invloggamma(y; branch=:lo)), y)
    end
end
