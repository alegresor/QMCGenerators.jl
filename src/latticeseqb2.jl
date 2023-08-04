mutable struct LatticeSeqB2
    const name::String
    const s::Int64 # dimension 
    const z::Vector{UInt64}
    const m::Int64 # 2^m points supported
    const n::Int64 # 2^m
    const scale::Float64 # 2^(-m)
    k::Int64 # index in the sequence
    v::UInt64 # current bitreversed binary
end

function LatticeSeqB2(s::Int64,z::Vector{BigInt},m::Int64)
    @assert m > 0 && s <= length(z)
    @assert m<=53
    n = 2^m 
    scale = Float64(2)^(-m)
    t = maximum(ndigits.(z[1:s],base=2))
    @assert t<=64
    z = convert.(UInt64,z[1:s])
    LatticeSeqB2("Lattice Seq B2",s,z,m,n,scale,-1,UInt64(0))
end

LatticeSeqB2(s::Int64,path::String,m::Int64) = LatticeSeqB2(s,readdlm(download(joinpath("https://bitbucket.org/dnuyens/qmc-generators/raw/cb0f2fb10fa9c9f2665e41419097781b611daa1e/LATSEQ/",path)),BigInt)[:,1],m)

LatticeSeqB2(s::Int64) = LatticeSeqB2(s,DEFAULT_LATTICESEQB2_GVECTOR,20)

function Reset!(seq::LatticeSeqB2)
    seq.k = -1
    seq.v = 0
    return
end

function NextLow(i1::UInt64,v::UInt64,xb::Matrix{Float64},s::Int64,n::Int64,k::Int64,z::Vector{UInt64},scale::Float64,m::Int64)
    b::Int64=0; frac::Float64=0
    for i::UInt64=i1:n
        k += 1
        b = 0; while ~Bool((k>>b)&1) b+= 1 end
        v ⊻= UInt64(1)<<(m-1-b)
        frac = scale*v
        for j=1:s xb[i,j] = convert(Float64,(frac*z[j])%1) end 
    end
    v
end 

function Next(seq::LatticeSeqB2,n::Int64)
    (seq.k+n)>=seq.n && throw(DomainError(n,"Generating $n more points will exceed the maximum supported points $(seq.n)"))
    xb::Matrix{Float64} = zeros(Float64,n,seq.s)
    seq.v = NextLow(seq.k==-1 ? UInt64(2) : UInt64(1),seq.v,xb,seq.s,n,seq.k==-1 ? 0 : seq.k,seq.z,seq.scale,seq.m)
    seq.k += n
    xb
end 

function FirstLinearLow(xb::Matrix{Float64},s::Int64,n::Int64,k::Int64,z::Vector{UInt64},scale::Float64,m::Int64)
    b::Int64=0; v::UInt64=0; frac::Float64=0
    for i::UInt64=1:n-1
        igc = (i⊻(i>>1))+1
        k += 1
        b = 0; while ~Bool((k>>b)&1) b+= 1 end
        v ⊻= UInt64(1)<<(m-1-b)
        frac = scale*v
        for j=1:s xb[igc,j] = convert(Float64,(frac*z[j])%1) end 
    end
end

function FirstLinear(seq::LatticeSeqB2,n::Int64)
    @assert seq.k == -1 
    @assert ispow2(n)
    (seq.k+n)>=seq.n && throw(DomainError(n,"Generating $n more points will exceed the maximum supported points $(seq.n)"))
    xb::Matrix{Float64} = zeros(Float64,n,seq.s)
    FirstLinearLow(xb,seq.s,n,seq.k==-1 ? 0 : seq.k,seq.z,seq.scale,seq.m)
    Reset!(seq)
    xb
end 

struct RandomShift
    name::String
    seq::LatticeSeqB2
    r::Int64
    rshifts::Matrix{Float64}
end 

RandomShift(seq::LatticeSeqB2,r::Int64,rng::Xoshiro) = RandomShift("Rand Shift: "*seq.name,seq,r,rand(rng,r,seq.s))

RandomShift(seq::LatticeSeqB2,r::Int64,seed::Int64) = RandomShift(seq,r,Xoshiro(seed))

RandomShift(seq::LatticeSeqB2,r::Int64) = RandomShift(seq,r,Xoshiro())

RandomShift(seq::LatticeSeqB2) = RandomShift(seq,1)

function NextRLowRandomShift(s::Int64,r::Int64,n::Int64,xu::Matrix{Float64},xr::Vector{Matrix{Float64}},rshifts::Matrix{Float64})
    for k=1:r 
        for j=1:s 
            for i=1:n
                xr[k][i,j] = (xu[i,j] + rshifts[k,j])%1
            end 
        end 
    end 
end

function NextR(rls::RandomShift,n::Int64)
    xu = Next(rls.seq,n)
    xr = [zeros(Float64,n,rls.seq.s) for k=1:rls.r]
    NextRLowRandomShift(rls.seq.s,rls.r,n,xu,xr,rls.rshifts)
    xr
end 

function FirstRLinear(rls::RandomShift,n::Int64)
    xu = FirstLinear(rls.seq,n)
    xr = [zeros(Float64,n,rls.seq.s) for k=1:rls.r]
    NextRLowRandomShift(rls.seq.s,rls.r,n,xu,xr,rls.rshifts)
    xr
end
