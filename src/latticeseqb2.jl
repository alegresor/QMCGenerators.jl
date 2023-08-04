mutable struct LatticeSeqB2
    const name::String
    const s::Int64 # dimension 
    const z::Vector{UInt64}
    const m::Int64 # 2^m points supported
    const n::Int64 # 2^m
    const scale::Union{BigFloat,Float64} # 2^(-m)
    k::Int64 # index in the sequence
end

function LatticeSeqB2(s::Int64,z::Vector{BigInt},m::Int64)
    @assert m > 0 && s <= length(z)
    n = 2^m 
    scale = m>53 ? BigFloat(2)^(-m) : Float64(2)^(-m)
    t = maximum(ndigits.(z[1:s],base=2))
    @assert t<=64
    z = convert.(UInt64,z[1:s])
    LatticeSeqB2("Lattice Seq B2",s,z,m,n,scale,-1)
end

LatticeSeqB2(s::Int64,path::String,m::Int64) = LatticeSeqB2(s,readdlm(download(joinpath("https://bitbucket.org/dnuyens/qmc-generators/raw/cb0f2fb10fa9c9f2665e41419097781b611daa1e/LATSEQ/",path)),BigInt)[:,1],m)

LatticeSeqB2(s::Int64) = LatticeSeqB2(s,DEFAULT_LATTICESEQB2_GVECTOR,20)

function Reset!(seq::LatticeSeqB2)
    seq.k = -1
    return
end

function NextLow(s::Int64,xb::Matrix{Float64},k::Int64,n::Int64,z::Vector{UInt64})
    idx::Int64,p2::Int64,frac::Float64 = 0,0,0
    for i=1:n 
        idx = k+i 
        p2 = 2^ceil(Int64,log2(idx+1))
        frac = (2*(p2-idx)-1)/p2
        for j=1:s
            xb[i,j] = convert(Float64,(frac*z[j])%1)
        end 
    end 
end 

function Next(seq::LatticeSeqB2,n::Int64)
    (seq.k+n)>=seq.n && throw(DomainError(n,"Generating $n more points will exceed the maximum supported points $(seq.n)"))
    xb::Matrix{Float64} = zeros(Float64,n,seq.s)
    NextLow(seq.s,xb,seq.k,n,seq.z)
    seq.k += n
    xb
end 

function FirstLinearLow(xb::Matrix{Float64},s::Int64,n::Int64,z::Vector{UInt64})
    frac::Float64 = 0
    for i=1:n
        frac = (i-1)/n
        for j=1:s
            xb[i,j] = convert(Float64,(frac*z[j])%1)
        end 
    end
end

function FirstLinear(seq::LatticeSeqB2,m::Int64)
    @assert seq.k == -1 
    n = 2^m
    (seq.k+n)>=seq.n && throw(DomainError(n,"Generating $n more points will exceed the maximum supported points $(seq.n)"))
    xb::Matrix{Float64} = zeros(Float64,n,seq.s)
    FirstLinearLow(xb,seq.s,n,seq.z)
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

function FirstRLinear(rls::RandomShift,m::Int64)
    n = 2^m
    xu = FirstLinear(rls.seq,m)
    xr = [zeros(Float64,n,rls.seq.s) for k=1:rls.r]
    NextRLowRandomShift(rls.seq.s,rls.r,n,xu,xr,rls.rshifts)
    xr
end