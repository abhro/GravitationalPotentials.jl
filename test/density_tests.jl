# SPDX-License-Identifier: MIT

using GravitationalPotentials
using Test
using Unitful, UnitfulAstro

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

        # TODO test model with scale radius == sphere radius
        r = 1u"pc"
        model = PowerLawSphereDensity(r, 1u"g/cm^3", r, -4)

        # TODO test model with scale radius < sphere radius
        model = PowerLawSphereDensity(r, 1u"g/cm^3", r/5, -4)

        # TODO test model with scale radius > sphere radius
        model = PowerLawSphereDensity(r, 1u"g/cm^3", r*5, -4)


        # TODO test model with α < -1
        model = PowerLawSphereDensity(1u"pc", 1u"g/cm^3", 0.5u"pc", -4)

        # TODO test model with α = -1
        model = PowerLawSphereDensity(1u"pc", 1u"g/cm^3", 0.5u"pc", -1)

        # TODO test model with α ∈ (-1, 0)
        model = PowerLawSphereDensity(1u"pc", 1u"g/cm^3", 0.5u"pc", -0.8)
        model = PowerLawSphereDensity(1u"pc", 1u"g/cm^3", 0.5u"pc", -0.5)

        # create model with α = 0
        @test_throws DimensionError PowerLawSphereDensity(1u"pc", 1u"g/cm^3", 0.5u"pc", 0)

        # TODO test model with α ∈ (0, 1)
        model = PowerLawSphereDensity(1u"pc", 1u"g/cm^3", 0.5u"pc", 0.5)
        model = PowerLawSphereDensity(1u"pc", 1u"g/cm^3", 0.5u"pc", 0.8)

        # TODO test model with α = 1
        model = PowerLawSphereDensity(1u"pc", 1u"g/cm^3", 0.5u"pc", 1)

        # TODO test model with α > 1
        model = PowerLawSphereDensity(1u"pc", 1u"g/cm^3", 0.5u"pc", 1.2)
        model = PowerLawSphereDensity(1u"pc", 1u"g/cm^3", 0.5u"pc", 4)
    end

    @testset "spiral galaxy" begin

    end

    @testset "Einasto galaxy" begin

    end

    @testset "NFW" begin
    end
end
