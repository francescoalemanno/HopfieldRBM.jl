using Documenter, HopfieldRBM

makedocs(
    modules = [HopfieldRBM],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Francesco Alemanno",
    sitename = "HopfieldRBM.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/francescoalemanno/HopfieldRBM.jl.git",
    push_preview = true
)
