include("../src/QMCGenerators.jl")
using .QMCGenerators
using PyPlot
using Random
import Distributions: quantile, Normal
import LinearAlgebra: norm

include("./plot.logo.jl")
plot_logo()

include("./plot.extensible_seq.jl")
plot_extensible_seq()

include("./plot.mc_vs_qmc.jl")
plot_mc_vs_qmc()
