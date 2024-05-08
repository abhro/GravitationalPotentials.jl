# SPDX-License-Identifier: MIT

using GravitationalPotentials
using Test

@testset "density models" begin

    @testset "uniform sphere" begin

        # TODO create model

        # TODO test at origin

        # TODO test at half radius

        # TODO test at random point between 0 and half radius

        # TODO test at random point between half and full radius

        # TODO test at radius/surface boundary

        # TODO test at 1.1 radius

        # TODO test at 2 radius

        # TODO test at infinity

    end

    @testset "uniform cylinder" begin

        # TODO create model

        # TODO test at roof

        # TODO test at floor

        # TODO test at walls

        # TODO test at origin

    end

    @testset "power law sphere" begin

        # TODO create model with scale radius == sphere radius
        # TODO create model with scale radius < sphere radius
        # TODO create model with scale radius > sphere radius


        # TODO create model with α < -1
        # TODO create model with α = -1
        # TODO create model with α ∈ (-1, 0)
        # TODO create model with α ∈ (0, 1)
        # TODO create model with α = 1
        # TODO create model with α > 1
    end

    @testset "spiral galaxy" begin

    end

    @testset "Einasto galaxy" begin

    end
end