push!(LOAD_PATH, "../src")

using Test
using Newcomb
using Printf


@testset "Newcomb" begin

   @test Newcomb.sunXYZ(2451545.5,2000.0) ≈ [0.1857261 , -0.8859448 , -0.3841002] atol = 1E-7

end

# ---------------------------------------------------------------------------------------------------------------------------------#
