mutable struct LatticeSeqB2
    const name::String
    const s::Int64 # dimension 
    const z::Union{Vector{UInt64},Vector{UInt128},Vector{BigInt}}
    const m::Int64 # 2^m points supported
    const n::Int64 # 2^m
    const scale::Union{BigFloat,Float64} # 2^(-m)
    const TInt::DataType
    k::Int64 # index in the sequence
end

function LatticeSeqB2(s::Int64,z::Vector{BigInt},m::Int64)
    @assert m > 0 && s <= length(z)
    n = 2^m 
    scale = m>53 ? BigFloat(2)^(-m) : Float64(2)^(-m)
    t = maximum(ndigits.(z[1:s],base=2))
    TInt = nothing; if t<=64 TInt = UInt64 elseif t<=128 TInt = UInt128 else TInt = BigInt end
    z = convert.(TInt,z[1:s])
    LatticeSeqB2("Lattice Seq B2",s,z,m,n,scale,TInt,-1)
end

LatticeSeqB2(s::Int64,path::String,m::Int64) = LatticeSeqB2(s,readdlm(download(joinpath("https://bitbucket.org/dnuyens/qmc-generators/raw/cb0f2fb10fa9c9f2665e41419097781b611daa1e/LATSEQ/",path)),BigInt)[:,1],m)

LatticeSeqB2(s::Int64) = LatticeSeqB2(s,DEFAULT_LATTICESEQB2_GVECTOR,20)

function Reset!(seq::LatticeSeqB2)
    seq.k = -1
    return
end

function NextLow(xb::Matrix{Float64},k::Int64,n::Int64,z::Union{Vector{UInt64},Vector{UInt128},Vector{BigInt}})
    nvec = range(k+1,k+n)
    p2s = 2 .^ ceil.(log2.(nvec.+1))
    fracs = (2 .* (p2s.-nvec) .- 1) ./ p2s
    xb[:,:] .= convert.(Float64,(fracs*z').%1)
end 

function Next(seq::LatticeSeqB2,n::Int64)
    (seq.k+n)>=seq.n && throw(DomainError(n,"Generating $n more points will exceed the maximum supported points $(seq.n)"))
    xb::Matrix{Float64} = zeros(Float64,n,seq.s)
    NextLow(xb,seq.k,n,seq.z)
    seq.k += n
    xb
end 

function FirstLinear(seq::LatticeSeqB2,m::Int64)
    @assert seq.k == -1 
    n = 2^m
    (seq.k+n)>=seq.n && throw(DomainError(n,"Generating $n more points will exceed the maximum supported points $(seq.n)"))
    convert.(Float64,(range(0,n-1)*(seq.z'./n)).%1)
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

function NextR(rls::RandomShift,n::Int64)
    xu = Next(rls.seq,n)
    [(xu.+rls.rshifts[i,:]').%1 for i=1:rls.r]
end 

function FirstRLinear(rls::RandomShift,m::Int64)
    xu = FirstLinear(rls.seq,m)
    [(xu.+rls.rshifts[i,:]').%1 for i=1:rls.r]
end
