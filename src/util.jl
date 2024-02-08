function jump!(rng::Xoshiro) # https://discourse.julialang.org/t/how-to-spawn-child-rngs-with-random/102122/3?u=alegresor
    JUMP = [0x180ec6d33cfd0aba,0xd5a61266f0c9392c,0xa9582618e03fc9aa,0x39abdc4529b1661c]
    s0,s1,s2,s3 = zero(UInt64),zero(UInt64),zero(UInt64),zero(UInt64)
    for j in JUMP
        for b in 0:63
            if (j & one(UInt64) << b) > 0
                s0 ⊻= rng.s0; s1 ⊻= rng.s1; s2 ⊻= rng.s2; s3 ⊻= rng.s3
            end
            rand(rng)
        end
    end
    rng.s0,rng.s1,rng.s2,rng.s3 = s0,s1,s2,s3
    return rng
end
spawn(rng::Xoshiro) = jump!(copy(rng))
function spawn(rng::Xoshiro,n::Int64)
    rngs = Vector{Xoshiro}(undef,n)
    for i=1:n rngs[i] = rng = spawn(rng) end 
    return rngs 
end

BinaryToFloat64(xb::UInt64,recipd::Union{Float64,BigFloat}) = convert.(Float64,recipd*xb)