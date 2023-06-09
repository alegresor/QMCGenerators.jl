using Test
using QMCGenerators
using Random

@testset "Documentation Plots" begin 
    include("./plot.logo.jl")
    plot_logo(false)

    include("./plot.extensible_seq.jl")
    plot_extensible_seq(2 .^ [2,3,4], false)

    include("./plot.mc_vs_qmc.jl")
    plot_mc_vs_qmc(3, 2, false)
    end
