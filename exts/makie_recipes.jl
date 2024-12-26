using Makie

@recipe(DensityHeatmap) do scene
    Attributes(
    )
end

convert_arguments()
function Makie.plot!(model::PowerLawCylinderDensity)
    heatmap()
end
