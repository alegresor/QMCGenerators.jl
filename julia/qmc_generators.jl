module QMCGenerators

import DelimitedFiles: readdlm
import Random: Xoshiro,MersenneTwister

bitreverse(v::BigInt, pad::Int64) = parse(BigInt,reverse(string(v,base=2,pad=pad)),base=2)


include("digitalseq_b2g.jl")
export DigitalSeqB2G,RandomDigitalShift,BinaryToFloat64,NextBinary,NextRBinary,FirstRLinearBinary,FirstLinearBinary

include("latticeseq_b2.jl")
export LatticeSeqB2,RandomShift

export Next,NextR,Reset!,FirstRLinear,FirstLinear

end 
