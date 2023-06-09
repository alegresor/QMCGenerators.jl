# https://documenter.juliadocs.org/stable/man/guide/ 

using Documenter, QMCGenerators

makedocs(
    sitename = "QMCGenerators.jl",
    modules = [QMCGenerators],
    format = Documenter.HTML(
        prettyurls = false, # get(ENV, "CI", nothing) == "true"
        sidebar_sitename = true
        ), 
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        # "API" => [
        #     "Lattice" => "latticeseqb2.md",
        #     "Digital Net" => "digitalseqb2g.md"],
    ]
)
