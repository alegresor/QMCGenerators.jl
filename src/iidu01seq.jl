mutable struct IIDU01Seq
    const name::String 
    const s::Int64 # dimension 
    const rngs::Vector{Xoshiro}
    const states::Matrix{UInt64}
    k::Int64 # index in the sequence
end

IIDU01Seq(s::Int64,rngs::Vector{Xoshiro}) = IIDU01Seq("IID Uniform01 Seq",s,rngs,vcat([[rngs[j].s0,rngs[j].s1,rngs[j].s2,rngs[j].s3] for j=1:s]'...),-1)

IIDU01Seq(s::Int64,rng::Xoshiro) = IIDU01Seq(s,spawn(rng,s))

IIDU01Seq(s::Int64,seed::Int64) = IIDU01Seq(s,Xoshiro(seed))

IIDU01Seq(s::Int64) = IIDU01Seq(s,Xoshiro())

function Reset!(seq::IIDU01Seq)
    for j=1:seq.s seq.rngs[j].s0,seq.rngs[j].s1,seq.rngs[j].s2,seq.rngs[j].s3 = seq.states[j,:] end 
    return
end

function Next(seq::IIDU01Seq,n::Int64)
    x = zeros(Float64,n,seq.s)
    for i=1:n
        seq.k = seq.k+1
        for j=1:seq.s x[i,j] = rand(seq.rngs[j]) end
    end 
    x 
end
