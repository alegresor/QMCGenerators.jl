mutable struct LatticeSeqB2
    name::String
    s::Int64 # dimension 
    z::Vector{BigInt} # generating vector 
    m::Int64 # 2^m points supported
    n::Int64 # 2^m
    scale::Union{BigFloat,Float64} # 2^(-m)
    k::Int64 # index in the sequence
end

function LatticeSeqB2(s::Int64,z::Vector{BigInt},m::Int64)
    @assert m > 0 && s <= length(z)
    n = 2^m 
    scale = m>53 ? BigFloat(2)^(-m) : Float64(2)^(-m)
    k = -1
    LatticeSeqB2("Lattice Seq B2",s,z,m,n,scale,k)
end

LatticeSeqB2(s::Int64,path::String,m::Int64) = LatticeSeqB2(s,readdlm(download(joinpath("https://bitbucket.org/dnuyens/qmc-generators/raw/cb0f2fb10fa9c9f2665e41419097781b611daa1e/LATSEQ/",path)),BigInt)[:,1],m)

LatticeSeqB2(s::Int64) = LatticeSeqB2(s,DEFAULT_LATTICESEQB2_GVECTOR,20)

function Reset!(seq::LatticeSeqB2)
    seq.k = -1 
    return
end

function Next(seq::LatticeSeqB2,n::Int64)
    x = zeros(Float64,n,seq.s)
    for i=1:n
        seq.k = seq.k+1 
        seq.k == seq.n && throw(DomainError(seq.k,"already generated maximum number of points"))
        phik = bitreverse(BigInt(seq.k),seq.m)*seq.scale 
        for j=1:seq.s
            xij = phik*seq.z[j]
            x[i,j] = convert(Float64,xij-floor(xij))
        end
    end 
    x 
end 

function FirstLinear(seq::LatticeSeqB2,m::Int64)
    n = 2^m
    x = [i for i=0:(n-1)]*seq.z[1:seq.s]'./n
    x = x.-floor.(x)
    convert.(Float64,x)
end 

mutable struct RandomShift
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
