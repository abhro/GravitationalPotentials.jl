# SPDX-License-Identifier: MIT

using GravitationalPotentials
using Documenter

DocMeta.setdocmeta!(
    GravitationalPotentials,
    :DocTestSetup,
    :(using GravitationalPotentials);
    recursive=true
)

makedocs(;
    modules=[GravitationalPotentials],
    authors="Abhro R. and contributors",
    sitename="GravitationalPotentials.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API reference" => "api.md",
    ],
)
