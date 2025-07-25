using AtmosphericModels
using Pkg
if ("TestEnv" ∈ keys(Pkg.project().dependencies))
    if ! ("Documenter" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
end
using Documenter

DocMeta.setdocmeta!(AtmosphericModels, :DocTestSetup, :(using AtmosphericModels); recursive=true)

makedocs(;
    modules=[AtmosphericModels],
    authors="Uwe Fechner <uwe.fechner.msc@gmail.com> and contributors",
    repo="https://github.com/OpenSourceAWE/KiteUtils.jl/blob/{commit}{path}#{line}",
    sitename="AtmosphericModels.jl",
    checkdocs=:none,
    format=Documenter.HTML(;
      repolink = "https://github.com/OpenSourceAWE/AtmosphericModels.jl",
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://OpenSourceAWE.github.io/AtmosphericModels.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API"  => "api.md",
        "Settings" => "settings.md",
        "3D Wind Fields" => "wind_field.md",
    ],
)

deploydocs(;
    repo="github.com/OpenSourceAWE/AtmosphericModels.jl",
    devbranch="main",
    push_preview=true,
)
