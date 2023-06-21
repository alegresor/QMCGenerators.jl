mutable struct IIDU01Seq
    name::String 
    s::Int64 # dimension 
    k::Int64 # index in the sequence
    x::Vector{Float64}
    rng::MersenneTwister
    internalseed::Vector{UInt32}
end

IIDU01Seq(s::Int64,rng::MersenneTwister) = IIDU01Seq("IID Uniform01 Seq",s,-1,zeros(Float64,s),rng,rng.seed)

IIDU01Seq(s::Int64,seed::Int64) = IIDU01Seq(s,MersenneTwister(seed))

IIDU01Seq(s::Int64) = IIDU01Seq(s,MersenneTwister())

function Reset!(seq::IIDU01Seq)
    seed!(seq.rng,seq.internalseed)
    return
end

function Next(seq::IIDU01Seq)
    seq.k = seq.k+1
    seq.x = rand(seq.rng,Float64,seq.s)
end

function Next(seq::IIDU01Seq,n::Int64)
    x = zeros(Float64,n,seq.s)
    for i=1:n
        x[i,:] = Next(seq)
    end 
    x 
end
