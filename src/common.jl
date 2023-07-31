Reset!(rseq::Union{RandomShift,RandomDigitalShift,RandomOwenScramble}) = Reset!(rseq.seq)

Next(seq::Union{IIDU01Seq,LatticeSeqB2,DigitalSeqB2G}) = Next(seq,1)[1,:]

NextR(rseq::Union{RandomShift,RandomDigitalShift,RandomOwenScramble}) = NextR(rseq,1)

function Next(rseq::Union{RandomShift,RandomDigitalShift,RandomOwenScramble},n)
    rseq.r!=1 && throw(DomainError(rseq.r,"Next requires 1 randomization"))
    NextR(rseq,n)[1]
end 

Next(rseq::Union{RandomShift,RandomDigitalShift,RandomOwenScramble}) = Next(rseq,1)
NextBinary(rseq::Union{RandomShift,RandomDigitalShift,RandomOwenScramble}) = NextBinary(rseq,1)

function FirstLinear(rseq::Union{RandomShift,RandomDigitalShift,RandomOwenScramble},n)
    rseq.r!=1 && throw(DomainError(rseq.r,"Next requires 1 randomization"))
    FirstRLinear(rseq,n)[1]
end

BinaryToFloat64(xb::Union{UInt64,Vector{UInt64},Matrix{UInt64}},seq::Union{DigitalSeqB2G,RandomDigitalShift,RandomOwenScramble}) = BinaryToFloat64.(xb,seq.recipd)

BinaryToFloat64(xbs::Vector{Matrix{UInt64}},rds::Union{RandomDigitalShift,RandomOwenScramble}) = [BinaryToFloat64(xb,rds) for xb in xbs]

NextRBinary(rds::Union{RandomDigitalShift,RandomOwenScramble}) = NextRBinary(rds,1)

function NextBinary(rds::Union{RandomDigitalShift,RandomOwenScramble},n::Int64)
    rds.r != 1 && throw(DomainError(rds.r,"Next requires 1 randomization"))
    NextRBinary(rds,n)[1]
end

NextBinary(rds::Union{RandomDigitalShift,RandomOwenScramble}) = NextBinary(rds,1)

NextR(rds::Union{RandomDigitalShift,RandomOwenScramble},n::Int64) = BinaryToFloat64(NextRBinary(rds,n),rds)

function FirstLinearBinary(rds::Union{RandomDigitalShift,RandomOwenScramble},m::Int64)
    rds.r != 1 && throw(DomainError(rds.r,"Next requires 1 randomization"))
    FirstRLinearBinary(rds,m)[1]
end

FirstRLinear(rds::Union{RandomDigitalShift,RandomOwenScramble},m::Int64) = BinaryToFloat64(FirstRLinearBinary(rds,m),rds)
