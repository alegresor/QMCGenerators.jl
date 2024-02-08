module QMCGenerators

import DelimitedFiles: readdlm
import Random: Xoshiro
import CairoMakie
using LaTeXStrings

include("util.jl")
export spawn

include("digitalseqb2g_default_gmatrix.jl")
include("digitalseqb2g.jl")
export LinearMatrixScramble,DigitalSeqB2G,RandomDigitalShift,RandomOwenScramble,BinaryToFloat64,NextRBinary,FirstRLinearBinary

include("latticeseqb2_default_gvector.jl")
include("latticeseqb2.jl")
export LatticeSeqB2,RandomShift

include("iidu01seq.jl")
export IIDU01Seq

include("common.jl")
export Reset!,Next,NextR,NextBinary,FirstLinear,BinaryToFloat64,NextBinary,FirstLinearBinary,FirstRLinear

include("plots.jl")
export qmcscatter!,JULIA4LOGOCOLORS

end 
