# SPDX-License-Identifier: MIT

using GravitationalPotentials
using Test
using Aqua

@testset "GravitationalPotentials.jl" begin
    include("density_tests.jl")

    @testset "Code quality (Aqua)" begin
        Aqua.test_all(GravitationalPotentials)
    end
end
