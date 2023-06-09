# https://documenter.juliadocs.org/stable/man/guide/ 

using Documenter, QMCGenerators

readme = read(joinpath(@__DIR__(),"../README.md"),String)
readme = replace(readme,"https://alegresor.github.io/QMCGenerators.jl/stable/tutorial"=>"@ref")
readme = replace(readme,"./docs/src"=>".")
write(joinpath(@__DIR__(),"./src/index.md"),readme)

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
