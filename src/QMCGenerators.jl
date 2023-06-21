module QMCGenerators

import DelimitedFiles: readdlm
import Random: MersenneTwister
import CairoMakie
using LaTeXStrings

bitreverse(v::BigInt, pad::Int64) = parse(BigInt,reverse(string(v,base=2,pad=pad)),base=2)

include("digitalseqb2g_default_gmatrix.jl")
include("digitalseqb2g.jl")
export DigitalSeqB2G,RandomDigitalShift,BinaryToFloat64,NextBinary,NextRBinary,FirstRLinearBinary,FirstLinearBinary

include("latticeseqb2_default_gvector.jl")
include("latticeseqb2.jl")
export LatticeSeqB2,RandomShift

include("util.jl")
export qmcscatter!,JULIA4LOGOCOLORS

export Next,NextR,Reset!,FirstRLinear,FirstLinear

end 
