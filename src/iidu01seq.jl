mutable struct IIDU01Seq
    name::String 
    s::Int64 # dimension 
    k::Int64 # index in the sequence
    x::Vector{Float64}
    rngs::Vector{Xoshiro}
    states::Matrix{UInt64}
end

IIDU01Seq(s::Int64,rngs::Vector{Xoshiro}) = IIDU01Seq("IID Uniform01 Seq",s,-1,zeros(Float64,s),rngs,vcat([[rngs[j].s0,rngs[j].s1,rngs[j].s2,rngs[j].s3] for j=1:s]'...))

IIDU01Seq(s::Int64,rng::Xoshiro) = IIDU01Seq(s,vcat(rng,[spawn(rng) for j=1:s-1]))

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
        seq.x = [rand(seq.rngs[j]) for j=1:seq.s]
        x[i,:] = seq.x
    end 
    x 
end
