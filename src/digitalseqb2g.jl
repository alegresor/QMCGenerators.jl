mutable struct DigitalSeqB2G
    name::String
    s::Int64 # dimension
    Csr::Matrix{BigInt} # matrix of size at least (s,m)
    m::Int64 # number of columns, can generate 2^m points
    t::Int64 # maximum number of bits in an element of Cs
    alpha::Float64 # t/m, the order of the net 
    n::Int64 # maximum number of supported points
    recipd::Union{BigFloat,Float64} # multiplication factor
    k::Int64 # index in the sequence
    cur::Array{BigInt}
end

function DigitalSeqB2G(s::Int64,Cs::Matrix{BigInt})
    @assert s <= size(Cs,1)
    m = size(Cs,2)
    t = maximum(ndigits.(Cs[1,:],base=2))
    Csr = bitreverse.(Cs[1:s,:],t)
    alpha = t/m 
    n = 2^m
    recipd = t>53 ? BigFloat(2)^(-t) : Float64(2)^(-t)
    cur = zeros(BigInt,s)
    DigitalSeqB2G("Digital Seq  B2",s,Csr,m,t,alpha,n,recipd,-1,cur)
end

DigitalSeqB2G(s::Int64,path::String) = DigitalSeqB2G(s,readdlm(download(joinpath("https://bitbucket.org/dnuyens/qmc-generators/raw/cb0f2fb10fa9c9f2665e41419097781b611daa1e/DIGSEQ/",path)),BigInt))

DigitalSeqB2G(s::Int64) = DigitalSeqB2G(s,DEFAULT_DIGITALSEQB2G_GMATRIX)

function Reset!(seq::DigitalSeqB2G) 
    seq.k = -1 
    seq.cur = zeros(BigInt,seq.s)
    return 
end 

function NextBinary!(seq::DigitalSeqB2G)
    seq.k += 1
    seq.k == seq.n && throw(DomainError(seq.k,"already generated maximum number of points"))
    if seq.k==0 return end 
    ctz = ndigits(((seq.k ⊻ (seq.k-1))+1) >> 1, base=2) - 1
    seq.cur .⊻= seq.Csr[:,(ctz+1)]
end

function NextBinary(seq::DigitalSeqB2G,n::Int64)
    xb = zeros(BigInt,n,seq.s)
    for i=1:n
        NextBinary!(seq)
        xb[i,:] = seq.cur
    end
    xb
end

BinaryToFloat64(xb::Union{BigInt,Vector{BigInt},Matrix{BigInt}},seq::DigitalSeqB2G) = convert.(Float64,seq.recipd*xb)

Next(seq::DigitalSeqB2G,n::Int64) = BinaryToFloat64(NextBinary(seq,n),seq)

Next(seq::DigitalSeqB2G) = Next(seq,1)

function FirstLinearBinary(seq::DigitalSeqB2G,m::Int64)
    @assert seq.k==-1
    n = 2^m
    xb = NextBinary(seq,n)
    gcs = map(i->i⊻(i>>1),0:n-1)
    xbl = zeros(BigInt,size(xb))
    for i=1:n xbl[gcs[i]+1,:] .= xb[i,:] end 
    Reset!(seq)
    xbl
end

FirstLinear(seq::DigitalSeqB2G,m::Int64) = BinaryToFloat64(FirstLinearBinary(seq,m),seq)

mutable struct RandomDigitalShift
    name::String
    seq::DigitalSeqB2G
    r::Int64
    rshifts::Matrix{BigInt}
    t::Int64 # number of bits in shifted integers
    tdiff::Int64
    recipd::Union{BigFloat,Float64} # multiplication factor
end 

function RandomDigitalShift(seq::DigitalSeqB2G,r::Int64,rng::MersenneTwister)
    t = max(seq.t,53)
    recipd = t>53 ? BigFloat(2)^(-t) : Float64(2)^(-t)
    rshifts = rand(rng,0:(BigInt(2)^t-1),r,seq.s)
    RandomDigitalShift("Digital Seq B2 + Random Shift",seq,r,rshifts,t,t-seq.t,recipd)
end

RandomDigitalShift(seq::DigitalSeqB2G,r::Int64,seed::Int64) = RandomDigitalShift(seq,r,MersenneTwister(seed))

RandomDigitalShift(seq::DigitalSeqB2G,r::Int64) = RandomDigitalShift(seq,r,MersenneTwister())

RandomDigitalShift(seq::DigitalSeqB2G) = RandomDigitalShift(seq,1)

function Reset!(rds::RandomDigitalShift) 
    Reset!(rds.seq)
end

DigitalShifts(xb,rds) = [rds.rshifts[i,:]' .⊻ xb for i=1:rds.r]

NextRBinary(rds::RandomDigitalShift,n::Int64) = DigitalShifts(NextBinary(rds.seq,n).<<rds.tdiff,rds)

NextRBinary(rds::RandomDigitalShift) = NextRBinary(rds,1)

function NextBinary(rds::RandomDigitalShift,n::Int64)
    rds.r != 1 && throw(DomainError(rds.r,"Next requires 1 randomization"))
    NextRBinary(rds,n)[1]
end

NextBinary(rds::RandomDigitalShift) = NextBinary(rds,1)

BinaryToFloat64(xb::Union{BigInt,Vector{BigInt},Matrix{BigInt}},rds::RandomDigitalShift) = convert.(Float64,rds.recipd*xb) # consolidate

BinaryToFloat64(xbs::Vector{Matrix{BigInt}},rds::RandomDigitalShift) = [BinaryToFloat64(xb,rds) for xb in xbs] # consolidate

NextR(rds::RandomDigitalShift,n::Int64) = BinaryToFloat64(NextRBinary(rds,n),rds)

NextR(rds::RandomDigitalShift) = NextR(rds,1)

function Next(rds::RandomDigitalShift,n)
    rds.r != 1 && throw(DomainError(rds.r,"Next requires 1 randomization"))
    NextR(rds,n)[1]
end

Next(rds::RandomDigitalShift) = Next(rds,1)

FirstRLinearBinary(rds::RandomDigitalShift,m::Int64) = DigitalShifts(FirstLinearBinary(rds.seq,m).<<rds.tdiff,rds)

function FirstLinearBinary(rds::RandomDigitalShift,m::Int64)
    rds.r != 1 && throw(DomainError(rds.r,"Next requires 1 randomization"))
    FirstRLinearBinary(rds,m)[1]
end 

FirstRLinear(rds::RandomDigitalShift,m::Int64) = BinaryToFloat64(FirstRLinearBinary(rds,m),rds)

function FirstLinear(rds::RandomDigitalShift,m::Int64)
    rds.r != 1 && throw(DomainError(rds.r,"Next requires 1 randomization"))
    FirstRLinear(rds,m)[1]
end
