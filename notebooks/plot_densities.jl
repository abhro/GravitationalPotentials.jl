### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 11757236-0cad-11ef-3235-69829e674f28
import Pkg; Pkg.activate(Base.current_project())

# ╔═╡ 90748e84-e2c3-43ff-b61b-900f1fed2a4f
using PlutoUI

# ╔═╡ 0dcb41a1-9cf3-4e4d-81bf-fdc860efa65c
using PlutoUI: TableOfContents, Slider, with_terminal

# ╔═╡ 8b711cb6-5a1f-4831-a136-ab6f1b96b1ec
using CairoMakie

# ╔═╡ 8e402635-93e8-4d67-b173-dc7532d2574c
using GravitationalPotentials

# ╔═╡ 0bdb3651-c955-4fb1-839f-8590bf2b2641
md"""
Use a shared environment for all notebooks demonstrating GravitationalPotentials.jl
"""

# ╔═╡ 301c7650-f3b9-4acf-9674-f4a4aa45e5d9
TableOfContents()

# ╔═╡ 59ed98e8-4d5b-490c-b0f6-5dd5ec363886
function density_heatmap(model, srange, zrange; title="Density heatmap", kwds...)
	fig = Figure()
	ax = Axis(fig[1,1], xlabel="s →", ylabel="z →"; title)

	domain = Iterators.product(srange, 0.0,  zrange)
	densities = mass_density.(model, domain)

	hm = heatmap!(ax, srange, zrange, densities; kwds...)
	cbar = Colorbar(fig[1,2], hm)
	return fig
end

# ╔═╡ 9f9d2922-8084-4a6c-9487-a10adb86373d
md"""
# Uniform density sphere
"""

# ╔═╡ 0780b58b-fa92-4c45-8b86-98c24683ec7b
@doc mass_density(::UniformSphereDensity, ::Any, ::Any, ::Any)

# ╔═╡ 0ceefa8f-e1dd-4b16-a22c-b5ce4e8f3f31
@bind sphere_radius NumberField(0:0.1:1000, default=0.5)

# ╔═╡ 7deda84f-1401-48db-980a-8c61e092b2bf
sphere_srange = range(0, 2.4, length=1000) * sphere_radius

# ╔═╡ 4af0e2fa-a05e-4350-bb86-71e5f8cdd4c7
sphere_zrange = range(-1.2, 1.2, length=100) * sphere_radius

# ╔═╡ 1ceee843-941f-4fd6-af44-24e891635a96
@bind sphere_density NumberField(0:0.1:1000, default=1.9)

# ╔═╡ bd760423-b837-4e3d-ab29-e0400baba63b
sphere_model = UniformSphereDensity(rₛ=sphere_radius, ρₛ=sphere_density)

# ╔═╡ 7f61c30b-145e-462e-ae01-8c05ebc440de
density_heatmap(sphere_model, sphere_srange, sphere_zrange)

# ╔═╡ 24ddf139-dec8-4c37-b29d-6fcbd6594662
md"""
# Power law sphere
"""

# ╔═╡ 1043fc7e-20dc-411c-80b8-3cfa909762c1
@doc mass_density(::PowerLawSphereDensity, ::Any, ::Any, ::Any)

# ╔═╡ df25a618-046e-48a1-918e-19d3978ef296
md"""
## Calculate density
"""

# ╔═╡ 213b4a3d-187e-4023-a4fb-fb0d4886f7f7
md"""
## Control parameters
Sphere radius: $(@bind pls_radius NumberField(0:0.1:1000, default=0.5)) $br
Scale radius: $(@bind pls_scale_radius NumberField(0:0.1:1000, default=0.5)) $br
Density: $(@bind pls_density NumberField(0:0.1:1000, default=1.9)) $br
α (exponent): $(@bind pls_alpha Slider(range(-5.01, 5.01, length=1000), default=1.5, show_value=true))
"""

# ╔═╡ b033332a-546f-4786-abd5-68300cebc1d8
pls_model = PowerLawSphereDensity{Float64,Float64}(
	rₛ = pls_radius,
	ρ₀ = pls_density,
	r₀ = pls_scale_radius,
	α = pls_alpha
)

# ╔═╡ 7d9ec66e-31ad-4f9f-8d24-a986b88dd8ba
pls_srange = range(0, 3, length=1000) * pls_radius

# ╔═╡ 6681ae1b-87e7-460b-9476-381345e49148
pls_zrange = range(-1.5, 1.5, length=100) * pls_radius

# ╔═╡ 33405c03-a162-4eac-9dd0-034628a0cbfe
md"""
## Plot power law sphere density model
"""

# ╔═╡ f1b4cf1a-5c92-4ae5-be96-0a4fd3374823
density_heatmap(pls_model, pls_srange, pls_zrange, title="Power law sphere density model", colormap=:jet)

# ╔═╡ 197d20f9-ca47-484c-b61b-654b82d25589
md"""
# Uniform density cylinder
"""

# ╔═╡ cba46d67-b50a-4afb-9b4c-76851501cd63
@doc mass_density(::UniformCylinderDensity, ::Any, ::Any, ::Any)

# ╔═╡ 0e6b434b-a45c-4229-924b-840aa20ef4da
md"""
## Calculate density
"""

# ╔═╡ 18f75cfc-6266-4fc7-9d67-af20d4ca9ab8
md"""
## Control cylinder model parameters
Radius: $(@bind cylinder_radius NumberField(0:0.1:1000, default=2)) $br
Height: $(@bind cylinder_height NumberField(0:0.1:1000, default=0.5)) $br
Density: $(@bind cylinder_density NumberField(0:0.1:1000, default=0.5)) $br
"""

# ╔═╡ c6c0f77f-0820-4516-90f1-4c9b6c182dca
cylinder_model = UniformCylinderDensity(
	r_c = cylinder_radius,
	h_c = cylinder_height,
	ρ_c = cylinder_density)

# ╔═╡ c9424081-3542-47d8-a4bc-a2945e264203
cyl_srange = range(0, 1.5, length=1000) * cylinder_radius

# ╔═╡ 29ad731e-68d7-4615-b329-7f299c637a7f
cyl_zrange = range(-1.5, 1.5, length=100) * cylinder_height

# ╔═╡ 9f3375fd-3f5f-470b-962a-4d27db8beb70
md"""
## Plot cylinder density model
"""

# ╔═╡ 6b9f6955-8637-4f9a-9fe5-064e5534dbed
density_heatmap(
    cylinder_model,
	cyl_srange, cyl_zrange,
    title="Uniform density cylinder model",
)

# ╔═╡ 2801a802-33cd-4617-983a-2c66d530ebf1
md"""
# Power law cylinder
"""

# ╔═╡ 8f01c85b-98d9-4e97-ba6a-8b09dbff1cfc
@doc mass_density(::PowerLawCylinderDensity, ::Any, ::Any, ::Any)

# ╔═╡ 0f1ec854-0011-4c3f-9a8b-ecf59faa87a7
# TODO

# ╔═╡ 3f975e70-3c53-4acc-9711-a5b47e31f1fe
md"""
# Spiral galaxy
"""

# ╔═╡ ba03a17e-5bdb-4138-aa23-b891724aa0db
@doc mass_density(::SpiralGalaxyDensity, ::Any, ::Any, ::Any)

# ╔═╡ fd2bd6b9-c1f8-469b-8fe0-aa4e7a12cec6
md"""
## Calculate galaxy density
"""

# ╔═╡ c2d4de7f-984f-48f8-ae33-789b185e4421
md"""
## Control galaxy shape
Bulge radius: $(@bind sg_bulge_radius NumberField(0:0.1:1000, default=1)) $br
Disk radius: $(@bind sg_disk_radius NumberField(0:0.1:1000, default=12)) $br
Disk height: $(@bind sg_disk_height NumberField(0:0.1:1000, default=0.5)) $br
"""

# ╔═╡ eaaec313-f8c8-42fc-90be-682248d5b6b9
sg_srange = range(0, 1.25, length=1000) * max(sg_bulge_radius, sg_disk_radius)

# ╔═╡ 63529ca5-4630-422a-92d2-60f86a8163ca
sg_zrange = range(-1.25, 1.25, length=1000) * max(sg_bulge_radius, sg_disk_height)

# ╔═╡ 50b1cf85-e197-4e69-b112-f3c88eeb35e9
md"""
## Control galaxy density
Bulge density: $(@bind sg_bulge_density NumberField(0:0.1:1000, default=190)) $br
Disk density: $(@bind sg_disk_density NumberField(0:0.1:1000, default=30)) $br
"""

# ╔═╡ ac6bbf27-68fc-49a1-8e91-28cedb283747
sg_model = SpiralGalaxyDensity(
    r_bulge=sg_bulge_radius,
    r_disk=sg_disk_radius,
    h_disk=sg_disk_height,

    ρ_bulge=sg_bulge_density,
    ρ_disk=sg_disk_density,
)

# ╔═╡ 06cc023c-0b5f-4755-a563-57011f198702
md"""
## Plot galaxy density
"""

# ╔═╡ d2dd3c5e-7feb-4068-857d-1ef31e51515b
density_heatmap(
	sg_model,
    sg_srange, sg_zrange,
    title="Spiral galaxy simplified density model",
    #colormap=:pastel,
)

# ╔═╡ 344682c5-5240-4e6d-bd66-6bfb6b58195a
md"""
# Einasto density model
"""

# ╔═╡ 75215ba3-2dd2-478a-9423-c19cbb7d53b9
@doc mass_density(::EinastoDensity, ::Any, ::Any, ::Any)

# ╔═╡ 0e7487a2-210f-43c3-82a5-53fdb96dc344
md"""
Central density (``ρ₀``): $(@bind einasto_ρ₀ NumberField(0:0.1:1000, default=1)) $br
Harmonic mean radius (``a₀``): $(@bind einasto_a₀ NumberField(0:0.1:1000, default=1)) $br
Structural parameter (``N``): $(@bind einasto_N NumberField(0:0.1:1000, default=1)) $br
Normalizing constant (``k``): $(@bind einasto_k NumberField(0:0.1:1000, default=1)) $br
"""

# ╔═╡ fe2d2ddc-64ff-42fc-9cc9-5da375b63a81
einasto_srange = range(0, 3.0, length=1000)

# ╔═╡ 6baa582f-0d4e-44e4-8ff7-e795a308276f
einasto_zrange = range(-3.0, 3.0, length=1000)

# ╔═╡ ac20970f-2543-4cab-933a-ef4420d8e91e
einasto_model = EinastoDensity(
    ρ₀=einasto_ρ₀,
	a₀=einasto_a₀,
	k=einasto_k,
	N=einasto_N,
)

# ╔═╡ 87b69f1b-5e68-45bb-bc6c-0db881fbcea9
density_heatmap(
	einasto_model,
	einasto_srange, einasto_zrange,
	#xlabel="s →", ylabel="z →",
	#transpose=true,
    title="Einasto density model",
)

# ╔═╡ 3119ef35-d927-4964-a4b7-57e52411bfff
md"""
# Navarro–Frenk–White profile
"""

# ╔═╡ 45836a41-acc6-4be7-98ed-07cce924b0c9
@doc mass_density(::NFWDensity, ::Any, ::Any, ::Any)

# ╔═╡ 43339592-2d27-4dee-875a-68e0388df9fd
md"""
# Build information
"""

# ╔═╡ 41aed17f-4dec-4902-9db3-4c02a9a196c5
versioninfo |> with_terminal

# ╔═╡ d130ad8f-b41f-4849-b73c-f6d3315f08aa
Pkg.status |> with_terminal

# ╔═╡ Cell order:
# ╟─0bdb3651-c955-4fb1-839f-8590bf2b2641
# ╠═11757236-0cad-11ef-3235-69829e674f28
# ╠═90748e84-e2c3-43ff-b61b-900f1fed2a4f
# ╠═0dcb41a1-9cf3-4e4d-81bf-fdc860efa65c
# ╠═8b711cb6-5a1f-4831-a136-ab6f1b96b1ec
# ╠═301c7650-f3b9-4acf-9674-f4a4aa45e5d9
# ╠═8e402635-93e8-4d67-b173-dc7532d2574c
# ╠═59ed98e8-4d5b-490c-b0f6-5dd5ec363886
# ╟─9f9d2922-8084-4a6c-9487-a10adb86373d
# ╠═0780b58b-fa92-4c45-8b86-98c24683ec7b
# ╠═bd760423-b837-4e3d-ab29-e0400baba63b
# ╠═7deda84f-1401-48db-980a-8c61e092b2bf
# ╠═4af0e2fa-a05e-4350-bb86-71e5f8cdd4c7
# ╠═0ceefa8f-e1dd-4b16-a22c-b5ce4e8f3f31
# ╠═1ceee843-941f-4fd6-af44-24e891635a96
# ╠═7f61c30b-145e-462e-ae01-8c05ebc440de
# ╟─24ddf139-dec8-4c37-b29d-6fcbd6594662
# ╠═1043fc7e-20dc-411c-80b8-3cfa909762c1
# ╟─df25a618-046e-48a1-918e-19d3978ef296
# ╠═b033332a-546f-4786-abd5-68300cebc1d8
# ╠═7d9ec66e-31ad-4f9f-8d24-a986b88dd8ba
# ╠═6681ae1b-87e7-460b-9476-381345e49148
# ╟─213b4a3d-187e-4023-a4fb-fb0d4886f7f7
# ╟─33405c03-a162-4eac-9dd0-034628a0cbfe
# ╠═f1b4cf1a-5c92-4ae5-be96-0a4fd3374823
# ╟─197d20f9-ca47-484c-b61b-654b82d25589
# ╠═cba46d67-b50a-4afb-9b4c-76851501cd63
# ╟─0e6b434b-a45c-4229-924b-840aa20ef4da
# ╠═c6c0f77f-0820-4516-90f1-4c9b6c182dca
# ╠═c9424081-3542-47d8-a4bc-a2945e264203
# ╠═29ad731e-68d7-4615-b329-7f299c637a7f
# ╟─18f75cfc-6266-4fc7-9d67-af20d4ca9ab8
# ╟─9f3375fd-3f5f-470b-962a-4d27db8beb70
# ╠═6b9f6955-8637-4f9a-9fe5-064e5534dbed
# ╟─2801a802-33cd-4617-983a-2c66d530ebf1
# ╠═8f01c85b-98d9-4e97-ba6a-8b09dbff1cfc
# ╠═0f1ec854-0011-4c3f-9a8b-ecf59faa87a7
# ╟─3f975e70-3c53-4acc-9711-a5b47e31f1fe
# ╠═ba03a17e-5bdb-4138-aa23-b891724aa0db
# ╟─fd2bd6b9-c1f8-469b-8fe0-aa4e7a12cec6
# ╠═ac6bbf27-68fc-49a1-8e91-28cedb283747
# ╠═eaaec313-f8c8-42fc-90be-682248d5b6b9
# ╠═63529ca5-4630-422a-92d2-60f86a8163ca
# ╟─c2d4de7f-984f-48f8-ae33-789b185e4421
# ╟─50b1cf85-e197-4e69-b112-f3c88eeb35e9
# ╟─06cc023c-0b5f-4755-a563-57011f198702
# ╠═d2dd3c5e-7feb-4068-857d-1ef31e51515b
# ╟─344682c5-5240-4e6d-bd66-6bfb6b58195a
# ╠═75215ba3-2dd2-478a-9423-c19cbb7d53b9
# ╟─0e7487a2-210f-43c3-82a5-53fdb96dc344
# ╠═fe2d2ddc-64ff-42fc-9cc9-5da375b63a81
# ╠═6baa582f-0d4e-44e4-8ff7-e795a308276f
# ╠═ac20970f-2543-4cab-933a-ef4420d8e91e
# ╠═87b69f1b-5e68-45bb-bc6c-0db881fbcea9
# ╟─3119ef35-d927-4964-a4b7-57e52411bfff
# ╠═45836a41-acc6-4be7-98ed-07cce924b0c9
# ╟─43339592-2d27-4dee-875a-68e0388df9fd
# ╠═41aed17f-4dec-4902-9db3-4c02a9a196c5
# ╠═d130ad8f-b41f-4849-b73c-f6d3315f08aa
