# import Pkg; Pkg.activate(@__DIR__)
let modpath = abspath(joinpath(dirname(@__DIR__), "src"))
    modpath in LOAD_PATH || push!(LOAD_PATH, modpath)
end
# import Pkg; Pkg.instantiate()
# push!(LOAD_PATH, "../../src/")
using Documenter, LocalPoly

makedocs(;
    sitename="LocalPoly.jl",
    modules=[LocalPoly],
    pages=[
        "Home" => "index.md",
        "API Reference" => "reference.md",
    ],
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        highlights=["matlab", "stata"]
    ),
)

deploydocs(
    repo="github.com/jbshannon/LocalPoly.jl.git",
    push_preview=true,
)
