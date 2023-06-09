mutable struct LatticeSeqB2
    s::Int64 # dimension 
    z::Vector{BigInt} # generating vector 
    m::Int64 # 2^m points supported
    n::Int64 # 2^m
    scale::Union{BigFloat,Float64} # 2^(-m)
    k::Int64 # index in the sequence
    x::Vector{BigFloat} 
end

function LatticeSeqB2(s::Int64,z::Vector{BigInt},m::Int64)
    @assert m>0
    n = 2^m 
    scale = m>53 ? BigFloat(2)^(-m) : Float64(2)^(-m)
    k = -1
    x = zeros(BigFloat,s)
    LatticeSeqB2(s,z,m,n,scale,k,x)
end

LatticeSeqB2(s::Int64,path::String,m::Int64) = LatticeSeqB2(s,readdlm(download(joinpath("https://raw.githubusercontent.com/alegresor/QMCGenerators.jl/main/src/LATSEQ/",path)),BigInt)[:,1],m)

LatticeSeqB2(s::Int64) = LatticeSeqB2(s,DEFAULT_LATTICESEQB2_GVECTOR,20)

function Reset!(ls::LatticeSeqB2)
    ls.k = -1 
end

function Next(ls::LatticeSeqB2)
    ls.k = ls.k+1 
    ls.k == ls.n && throw(DomainError(ls.k,"already generated maximum number of points"))
    phik = bitreverse(BigInt(ls.k),ls.m)*ls.scale 
    for j=1:ls.s 
        ls.x[j] = phik*ls.z[j]
        ls.x[j] = ls.x[j]-floor(ls.x[j])
    end
    convert.(Float64,ls.x)
end

function Next(ls::LatticeSeqB2,n::Int64)
    x = zeros(Float64,n,ls.s)
    for i=1:n
        x[i,:] = Next(ls)
    end 
    x 
end 

function FirstLinear(ls::LatticeSeqB2,m::Int64)
    n = 2^m
    x = [i for i=0:(n-1)]*ls.z[1:ls.s]'./n
    x = x.-floor.(x)
    convert.(Float64,x)
end 

mutable struct RandomShift
    ls::LatticeSeqB2
    r::Int64
    rshifts::Matrix{Float64}
end 

RandomShift(ls::LatticeSeqB2,r::Int64,rng::MersenneTwister) = RandomShift(ls,r,rand(rng,r,ls.s))

RandomShift(ls::LatticeSeqB2,r::Int64,seed::Int64) = RandomShift(ls,r,MersenneTwister(seed))

RandomShift(ls::LatticeSeqB2,r::Int64) = RandomShift(ls,r,MersenneTwister())

RandomShift(ls::LatticeSeqB2) = RandomShift(ls,1)

function Reset!(rls::RandomShift)
    Reset!(rls.ls)
end

function NextR(rls::RandomShift,n::Int64)
    xu = Next(rls.ls,n)
    [(xu.+rls.rshifts[i,:]').%1 for i=1:rls.r]
end 

NextR(rls::RandomShift) = NextR(rls,1)

function Next(rls::RandomShift,n)
    rls.r!=1 && throw(DomainError(rls.r,"Next requires 1 randomization"))
    NextR(rls,n)[1]
end 

Next(rls::RandomShift) = Next(rls,1)

function FirstRLinear(rls::RandomShift,m::Int64)
    xu = FirstLinear(rls.ls,m)
    [(xu.+rls.rshifts[i,:]').%1 for i=1:rls.r]
end 

function FirstLinear(rls::RandomShift,n)
    rls.r!=1 && throw(DomainError(rls.r,"Next requires 1 randomization"))
    FirstRLinear(rls,n)[1]
end
