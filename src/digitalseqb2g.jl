mutable struct DigitalSeqB2G
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
    m = size(Cs,2)
    t = maximum(ndigits.(Cs[1,:],base=2))
    Csr = bitreverse.(Cs[1:s,:],t)
    alpha = t/m 
    n = 2^m
    recipd = t>53 ? BigFloat(2)^(-t) : Float64(2)^(-t)
    cur = zeros(BigInt,s)
    DigitalSeqB2G(s,Csr,m,t,alpha,n,recipd,-1,cur)
end

DigitalSeqB2G(s::Int64,path::String) = DigitalSeqB2G(s,readdlm(download(joinpath("https://raw.githubusercontent.com/alegresor/QMCGenerators.jl/main/src/DIGSEQ/",path)),BigInt))

DigitalSeqB2G(s::Int64) = DigitalSeqB2G(s,DEFAULT_DIGITALSEQB2G_GMATRIX)

function Reset!(ds::DigitalSeqB2G) 
    ds.k = -1 
    ds.cur = zeros(BigInt,ds.s)
end 

function NextBinary!(ds::DigitalSeqB2G)
    ds.k += 1
    ds.k == ds.n && throw(DomainError(ds.k,"already generated maximum number of points"))
    if ds.k==0 return end 
    ctz = ndigits(((ds.k ⊻ (ds.k-1))+1) >> 1, base=2) - 1
    ds.cur .⊻= ds.Csr[:,(ctz+1)]
end

function NextBinary(ds::DigitalSeqB2G,n::Int64)
    xb = zeros(BigInt,n,ds.s)
    for i=1:n
        NextBinary!(ds)
        xb[i,:] = ds.cur
    end
    xb
end

BinaryToFloat64(xb::Union{BigInt,Vector{BigInt},Matrix{BigInt}},ds::DigitalSeqB2G) = convert.(Float64,ds.recipd*xb)

Next(ds::DigitalSeqB2G,n::Int64) = BinaryToFloat64(NextBinary(ds,n),ds)

Next(ds::DigitalSeqB2G) = Next(ds,1)

function FirstLinearBinary(ds::DigitalSeqB2G,m::Int64)
    @assert ds.k==-1
    n = 2^m
    xb = NextBinary(ds,n)
    gcs = map(i->i⊻(i>>1),0:n-1)
    xbl = zeros(BigInt,size(xb))
    for i=1:n xbl[gcs[i]+1,:] .= xb[i,:] end 
    Reset!(ds)
    xbl
end

FirstLinear(ds::DigitalSeqB2G,m::Int64) = BinaryToFloat64(FirstLinearBinary(ds,m),ds)

mutable struct RandomDigitalShift
    ds::DigitalSeqB2G
    r::Int64
    rshifts::Matrix{BigInt}
end 

RandomDigitalShift(ds::DigitalSeqB2G,r::Int64,rng::MersenneTwister) = RandomDigitalShift(ds,r,rand(rng,0:(BigInt(2)^ds.t-1),r,ds.s))

RandomDigitalShift(ds::DigitalSeqB2G,r::Int64,seed::Int64) = RandomDigitalShift(ds,r,MersenneTwister(seed))

RandomDigitalShift(ds::DigitalSeqB2G,r::Int64) = RandomDigitalShift(ds,r,MersenneTwister())

RandomDigitalShift(ds::DigitalSeqB2G) = RandomDigitalShift(ds,1)

function Reset!(rds::RandomDigitalShift) 
    Reset!(rds.ds)
end

DigitalShifts(xb,rds) = [rds.rshifts[i,:]' .⊻ xb for i=1:rds.r]

NextRBinary(rds::RandomDigitalShift,n::Int64) = DigitalShifts(NextBinary(rds.ds,n),rds)

NextRBinary(rds::RandomDigitalShift) = NextRBinary(rds,1)

function NextBinary(rds::RandomDigitalShift,n::Int64)
    rds.r != 1 && throw(DomainError(rds.r,"Next requires 1 randomization"))
    NextRBinary(rds,n)[1]
end

NextBinary(rds::RandomDigitalShift) = NextBinary(rds,1)

BinaryToFloat64(xb::Union{BigInt,Vector{BigInt},Matrix{BigInt}},rds::RandomDigitalShift) = BinaryToFloat64(xb,rds.ds)

BinaryToFloat64(xbs::Vector{Matrix{BigInt}},rds::RandomDigitalShift) = [BinaryToFloat64(xb,rds) for xb in xbs]

NextR(rds::RandomDigitalShift,n::Int64) = BinaryToFloat64(NextRBinary(rds,n),rds)

NextR(rds::RandomDigitalShift) = NextR(rds,1)

function Next(rds::RandomDigitalShift,n)
    rds.r != 1 && throw(DomainError(rds.r,"Next requires 1 randomization"))
    NextR(rds,n)[1]
end

Next(rds::RandomDigitalShift) = Next(rds,1)

FirstRLinearBinary(rds::RandomDigitalShift,m::Int64) = DigitalShifts(FirstLinearBinary(rds.ds,m),rds)

function FirstLinearBinary(rds::RandomDigitalShift,m::Int64)
    rds.r != 1 && throw(DomainError(rds.r,"Next requires 1 randomization"))
    FirstRLinearBinary(rds,m)[1]
end 

FirstRLinear(rds::RandomDigitalShift,m::Int64) = BinaryToFloat64(FirstRLinearBinary(rds,m),rds)

function FirstLinear(rds::RandomDigitalShift,m::Int64)
    rds.r != 1 && throw(DomainError(rds.r,"Next requires 1 randomization"))
    FirstRLinear(rds,m)[1]
end 
