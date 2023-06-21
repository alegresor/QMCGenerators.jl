# https://documenter.juliadocs.org/stable/man/guide/ 

using Documenter, QMCGenerators

readme = read(joinpath(@__DIR__(),"../README.md"),String)
readme = replace(readme,"https://alegresor.github.io/QMCGenerators.jl/stable/tutorial"=>"@ref")
readme = replace(readme,"./docs/src"=>".")
write(joinpath(@__DIR__(),"./src/index.md"),readme)

assetsdir = joinpath(@__DIR__,"src/assets")
if ~isdir(assetsdir) mkdir(assetsdir) end

println("DOCTEST")
makedocs(
    doctest = :only,
    strict = true
)

println("BUILDING DOCS")
makedocs(
    sitename = "QMCGenerators.jl",
    modules = [QMCGenerators],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        sidebar_sitename = true
        ), 
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        # "API" => [
        #     "Lattice" => "latticeseqb2.md",
        #     "Digital Net" => "digitalseqb2g.md"],
    ],
    doctest = false
)

println("DEPLOY DOCS")
deploydocs(
    repo = "github.com/alegresor/QMCGenerators.jl.git",
    devbranch = "main",
    #versions = ["stable" => "v^", "main" => "main"]
)
