module QMCGenerators

import DelimitedFiles: readdlm
import Random: Xoshiro
import CairoMakie
using LaTeXStrings

include("util.jl")

include("digitalseqb2g_default_gmatrix.jl")
include("digitalseqb2g.jl")
export DigitalSeqB2G,RandomOwenScramble,LinearMatrixScramble,RandomDigitalShift,BinaryToFloat64,NextBinary,NextRBinary,FirstRLinearBinary,FirstLinearBinary

include("latticeseqb2_default_gvector.jl")
include("latticeseqb2.jl")
export LatticeSeqB2,RandomShift

include("iidu01seq.jl")
export IIDU01Seq

include("common.jl")
export Next,NextR,Reset!,FirstRLinear,FirstLinear

include("plots.jl")
export qmcscatter!,JULIA4LOGOCOLORS

end 
