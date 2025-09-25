using InverseGammaFunction
using Test

@testset "InverseGammaFunction.jl" begin
    @test isapprox(invgamma(100.0), 5.892518696343773)
    @test isapprox(invgamma(100.0; branch=:lo), 0.00994357323626315)
end
