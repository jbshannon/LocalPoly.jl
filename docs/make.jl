# import Pkg; Pkg.activate(@__DIR__)
let modpath = abspath(joinpath(dirname(@__DIR__), "src"))
    modpath in LOAD_PATH || push!(LOAD_PATH, modpath)
end
# import Pkg; Pkg.instantiate()
# push!(LOAD_PATH, "../../src/")
using Documenter, LocalPoly

makedocs(;
    sitename="LocalPoly.jl",
    pages=[
        "Home" => "index.md",
        "API Reference" => "reference.md",
    ]
)

deploydocs(
    repo="github.com/jbshannon/LocalPoly.jl.git",
    versions = ["stable" => "v^", "v#.#", "dev" => "master"],
    push_preview=true,
)
