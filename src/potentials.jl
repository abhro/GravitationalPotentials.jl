using Integrals
using LinearAlgebra: norm
function potential(model::MassDensityModel, r)
    # get the domain of the model
    domain = bounds_t(model) # FIXME we're using Extents.jl now

    # integrate over the model
    function f!(y, r′, p)
        y .= mass_density(model, r′) / norm(r′ - p)
    end
    prob = IntegralProblem(f, domain, r)
    return solve(prob, HCubatureJL())
end
