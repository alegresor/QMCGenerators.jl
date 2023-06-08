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

LatticeSeqB2(s::Int64,path::String,m::Int64) = LatticeSeqB2(s,readdlm(joinpath(@__DIR__(),"LATSEQ",path),BigInt)[:,1],m)

LatticeSeqB2(s::Int64) = LatticeSeqB2(s,
    BigInt[1, 182667, 469891, 498753, 110745, 446247, 250185, 118627, 245333, 283199, 
    408519, 391023, 246327, 126539, 399185, 461527, 300343, 69681, 516695, 436179, 106383, 238523, 
    413283, 70841, 47719, 300129, 113029, 123925, 410745, 211325, 17489, 511893, 40767, 186077, 
    519471, 255369, 101819, 243573, 66189, 152143, 503455, 113217, 132603, 463967, 297717, 157383, 
    224015, 502917, 36237, 94049, 170665, 79397, 123963, 223451, 323871, 303633, 98567, 318855, 
    494245, 477137, 177975, 64483, 26695, 88779, 94497, 239429, 381007, 110205, 339157, 73397, 
    407559, 181791, 442675, 301397, 32569, 147737, 189949, 138655, 350241, 63371, 511925, 515861, 
    434045, 383435, 249187, 492723, 479195, 84589, 99703, 239831, 269423, 182241, 61063, 130789, 
    143095, 471209, 139019, 172565, 487045, 304803, 45669, 380427, 19547, 425593, 337729, 237863, 
    428453, 291699, 238587, 110653, 196113, 465711, 141583, 224183, 266671, 169063, 317617, 68143,
    291637, 263355, 427191, 200211, 365773, 254701, 368663, 248047, 209221, 279201, 323179, 80217, 
    122791, 316633, 118515, 14253, 129509, 410941, 402601, 511437, 10469, 366469, 463959, 442841, 
    54641, 44167, 19703, 209585, 69037, 33317, 433373, 55879, 245295, 10905, 468881, 128617, 417919, 
    45067, 442243, 359529, 51109, 290275, 168691, 212061, 217775, 405485, 313395, 256763, 152537, 326437, 
    332981, 406755, 423147, 412621, 362019, 279679, 169189, 107405, 251851, 5413, 316095, 247945, 422489, 
    2555, 282267, 121027, 369319, 204587, 445191, 337315, 322505, 388411, 102961, 506099, 399801, 254381, 
    452545, 309001, 147013, 507865, 32283, 320511, 264647, 417965, 227069, 341461, 466581, 386241, 
    494585, 201479, 151243, 481337, 68195, 75401, 58359, 448107, 459499, 9873, 365117, 350845, 181873, 
    7917, 436695, 43899, 348367, 423927, 437399, 385089, 21693, 268793, 49257, 250211, 125071, 341631, 
    310163, 94631, 108795, 21175, 142847, 383599, 71105, 65989, 446433, 177457, 107311, 295679, 442763, 
    40729, 322721, 420175, 430359, 480757],20)

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

RandomShift(ls::LatticeSeqB2,r::Int64,rng::Union{Xoshiro,MersenneTwister}) = RandomShift(ls,r,rand(rng,r,ls.s))

RandomShift(ls::LatticeSeqB2,r::Int64,seed::Int64) = RandomShift(ls,r,Xoshiro(seed))

RandomShift(ls::LatticeSeqB2,r::Int64) = RandomShift(ls,r,Xoshiro())

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
