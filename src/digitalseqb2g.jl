struct LinearMatrixScramble
    s::Int64 
    Csrlms::Union{Matrix{UInt64},Matrix{UInt128},Matrix{BigInt}}
    m::Int64
    t::Int64
    TInt::DataType
end

function LinearMatrixScramble(s::Int64,Cs::Matrix{BigInt},rng::Xoshiro) # Matousek
    @assert s <= size(Cs,1)
    m = size(Cs,2)
    tog = maximum(ndigits.(Cs[1:s,:],base=2))
    t = max(53,tog)
    TInt = nothing; if t<=64 TInt = UInt64 elseif t<=128 TInt = UInt128 else TInt = BigInt end
    Csr = convert.(TInt,bitreverse.(Cs[1:s,:],tog))
    Csrlms = zeros(TInt,s,m)
    for j=1:s
        for k=0:t-1
            k1 = min(k,tog)
            delem = TInt(1)<<(tog-k1-1)
            ndelem = rand(rng,0:(TInt(1)<<k1)-1)<<(tog-k1)
            u = delem+ndelem
            for l=1:m
                v = u&Csr[j,l]
                setbits = 0 
                for i=0:tog setbits += (v>>i)&1 end 
                b = setbits%2
                Csrlms[j,l] += TInt(b)<<(t-k-1)
            end
        end
    end
    LinearMatrixScramble(s,Csrlms,m,t,TInt)
end

LinearMatrixScramble(s::Int64,Cs::Matrix{BigInt},seed::Int64) = LinearMatrixScramble(s,Cs,Xoshiro(seed))

LinearMatrixScramble(s::Int64,Cs::Matrix{BigInt}) = LinearMatrixScramble(s,Cs,Xoshiro())

LinearMatrixScramble(s::Int64,path::String,rng::Xoshiro) = LinearMatrixScramble(s,readdlm(download(joinpath("https://bitbucket.org/dnuyens/qmc-generators/raw/cb0f2fb10fa9c9f2665e41419097781b611daa1e/DIGSEQ/",path)),BigInt),rng)
LinearMatrixScramble(s::Int64,path::String,seed::Int64) = LinearMatrixScramble(s,readdlm(download(joinpath("https://bitbucket.org/dnuyens/qmc-generators/raw/cb0f2fb10fa9c9f2665e41419097781b611daa1e/DIGSEQ/",path)),BigInt),seed)
LinearMatrixScramble(s::Int64,path::String) = LinearMatrixScramble(s,readdlm(download(joinpath("https://bitbucket.org/dnuyens/qmc-generators/raw/cb0f2fb10fa9c9f2665e41419097781b611daa1e/DIGSEQ/",path)),BigInt))

LinearMatrixScramble(s::Int64,rng::Xoshiro) = LinearMatrixScramble(s,DEFAULT_DIGITALSEQB2G_GMATRIX,rng)
LinearMatrixScramble(s::Int64,seed::Int64) = LinearMatrixScramble(s,DEFAULT_DIGITALSEQB2G_GMATRIX,seed)
LinearMatrixScramble(s::Int64) = LinearMatrixScramble(s,DEFAULT_DIGITALSEQB2G_GMATRIX)

mutable struct DigitalSeqB2G
    const name::String
    const s::Int64 # dimension
    const Csr::Union{Matrix{UInt64},Matrix{UInt128},Matrix{BigInt}}
    const m::Int64 # number of columns, can generate 2^m points
    const t::Int64 # maximum number of bits in an element of Cs
    const alpha::Float64 # t/m, the order of the net 
    const n::Int64 # maximum number of supported points
    const recipd::Union{BigFloat,Float64} # multiplication factor
    const TInt::DataType
    k::Int64 # index in the sequence
    cur::Union{Vector{UInt64},Vector{UInt128},Vector{BigInt}}
end

function DigitalSeqB2G(s::Int64,Cs::Matrix{BigInt})
    @assert s <= size(Cs,1)
    m = size(Cs,2)
    t = maximum(ndigits.(Cs[1:s,:],base=2))
    alpha = t/m 
    n = 2^m
    recipd = t>53 ? BigFloat(2)^(-t) : Float64(2)^(-t)
    TInt = nothing; if t<=64 TInt = UInt64 elseif t<=128 TInt = UInt128 else TInt = BigInt end
    Csr = convert.(TInt,bitreverse.(Cs[1:s,:],t))
    cur = zeros(TInt,s)
    DigitalSeqB2G("Digital Seq B2",s,Csr,m,t,alpha,n,recipd,TInt,-1,cur)
end

DigitalSeqB2G(s::Int64,path::String) = DigitalSeqB2G(s,readdlm(download(joinpath("https://bitbucket.org/dnuyens/qmc-generators/raw/cb0f2fb10fa9c9f2665e41419097781b611daa1e/DIGSEQ/",path)),BigInt))

DigitalSeqB2G(s::Int64) = DigitalSeqB2G(s,DEFAULT_DIGITALSEQB2G_GMATRIX)

function DigitalSeqB2G(rlms::LinearMatrixScramble)
    alpha = rlms.t/rlms.m
    n = 2^rlms.m
    recipd = rlms.t>53 ? BigFloat(2)^(-rlms.t) : Float64(2)^(-rlms.t)
    cur = zeros(rlms.TInt,rlms.s)
    DigitalSeqB2G("LMS Digital Seq B2",rlms.s,rlms.Csrlms,rlms.m,rlms.t,alpha,n,recipd,rlms.TInt,-1,cur)
end 

function Reset!(seq::DigitalSeqB2G) 
    seq.k = -1 
    seq.cur = zeros(seq.TInt,seq.s)
    return 
end 

function NextBinaryLow(i1::Int64,s::Int64,xb::Union{Matrix{UInt64},Matrix{UInt128},Matrix{BigInt}},n::Int64,k::Int64,Csr::Union{Matrix{UInt64},Matrix{UInt128},Matrix{BigInt}},cur::Union{Vector{UInt64},Vector{UInt128},Vector{BigInt}})
    b::Int64 = 0
    for i::Int64=i1:n
        k += 1
        b = 0; while ~Bool((k>>b)&1) b+= 1 end
        for j=1:s xb[i,j] = cur[j] ⊻= Csr[j,b+1] end
    end
end

function NextBinary(seq::DigitalSeqB2G,n::Int64)
    (seq.k+n)>=seq.n && throw(DomainError(n,"Generating $n more points will exceed the maximum supported points $(seq.n)"))
    xb::Union{Matrix{UInt64},Matrix{UInt128},Matrix{BigInt}} = zeros(seq.TInt,n,seq.s)    
    NextBinaryLow(seq.k==-1 ? 2 : 1,seq.s,xb,n,seq.k==-1 ? 0 : seq.k,seq.Csr,seq.cur)
    seq.k += n 
    seq.cur .= xb[n,:]
    xb
end

Next(seq::DigitalSeqB2G,n::Int64) = BinaryToFloat64(NextBinary(seq,n),seq)

function FirstLinearBinaryLow(s::Int64,xb::Union{Matrix{UInt64},Matrix{UInt128},Matrix{BigInt}},n::Int64,k::Int64,Csr::Union{Matrix{UInt64},Matrix{UInt128},Matrix{BigInt}},cur::Union{Vector{UInt64},Vector{UInt128},Vector{BigInt}})
    b::Int64 = 0
    for i::Int64=1:n-1
        igc = (i⊻(i>>1))+1
        k += 1
        b = 0; while ~Bool((k>>b)&1) b+= 1 end
        for j=1:s xb[igc,j] = cur[j] ⊻= Csr[j,b+1] end
    end
end

function FirstLinearBinary(seq::DigitalSeqB2G,m::Int64)
    n = 2^m
    n>=seq.n && throw(DomainError(n,"Generating $n more points will exceed the maximum supported points $(seq.n)"))
    @assert seq.k==-1
    xb::Union{Matrix{UInt64},Matrix{UInt128},Matrix{BigInt}} = zeros(seq.TInt,n,seq.s)
    FirstLinearBinaryLow(seq.s,xb,n,seq.k==-1 ? 0 : seq.k,seq.Csr,seq.cur)
    Reset!(seq)
    xb
end

FirstLinear(seq::DigitalSeqB2G,m::Int64) = BinaryToFloat64(FirstLinearBinary(seq,m),seq)

struct RandomDigitalShift
    name::String
    seq::DigitalSeqB2G
    r::Int64
    rshifts::Union{Matrix{UInt64},Matrix{UInt128},Matrix{BigInt}}
    t::Int64 # number of bits in shifted integers
    tdiff::Int64
    recipd::Union{BigFloat,Float64} # multiplication factor
end 

function RandomDigitalShift(seq::DigitalSeqB2G,r::Int64,rng::Xoshiro)
    t = max(53,seq.t)
    recipd = t>53 ? BigFloat(2)^(-t) : Float64(2)^(-t)
    rshifts = rand(rng,0:(seq.TInt(2)^t-1),r,seq.s)
    RandomDigitalShift("Rand DShift: "*seq.name,seq,r,rshifts,t,t-seq.t,recipd)
end

RandomDigitalShift(seq::DigitalSeqB2G,r::Int64,seed::Int64) = RandomDigitalShift(seq,r,Xoshiro(seed))

RandomDigitalShift(seq::DigitalSeqB2G,r::Int64) = RandomDigitalShift(seq,r,Xoshiro())

RandomDigitalShift(seq::DigitalSeqB2G) = RandomDigitalShift(seq,1)

DigitalShifts(xb::Union{Matrix{UInt64},Matrix{UInt128},Matrix{BigInt}},rds::RandomDigitalShift) = [rds.rshifts[i,:]' .⊻ xb for i=1:rds.r]

NextRBinary(rds::RandomDigitalShift,n::Int64) = DigitalShifts(NextBinary(rds.seq,n).<<rds.tdiff,rds)

FirstRLinearBinary(rds::RandomDigitalShift,m::Int64) = DigitalShifts(FirstLinearBinary(rds.seq,m).<<rds.tdiff,rds)

mutable struct BTree
    rbits::Union{Nothing,UInt64,UInt128,BigInt,Bool}
    xb::Union{Nothing,UInt64,UInt128,BigInt}
    left::Union{Nothing,BTree}
    right::Union{Nothing,BTree}
end
BTree(left::BTree,right::BTree) = BTree(nothing,nothing,left,right)
BTree(rbits::Union{UInt64,UInt128,BigInt},xb::Union{UInt64,UInt128,BigInt}) = BTree(rbits,xb,nothing,nothing)
BTree() = BTree(nothing,nothing,nothing,nothing)

struct RandomOwenScramble
    name::String
    seq::DigitalSeqB2G
    r::Int64
    rngs::Matrix{Xoshiro}
    scrambles::Matrix{BTree}
    t::Int64 # number of bits in shifted integers
    tdiff::Int64
    recipd::Union{BigFloat,Float64} # multiplication factor
end

function RandomOwenScramble(seq::DigitalSeqB2G,r::Int64,rngs::Matrix{Xoshiro})
    t = max(53,seq.t)
    recipd = t>53 ? BigFloat(2)^(-t) : Float64(2)^(-t)
    scrambles = Matrix{BTree}(undef,r,seq.s)
    for i=1:r 
        for j=1:seq.s 
            r1 = seq.TInt(rand(rngs[i,j],Bool))<<(t-1)
            rbitsleft,rbitsright = r1 + rand(rngs[i,j],0:(seq.TInt(2)^(t-1)-1)),r1 + rand(rngs[i,j],0:(seq.TInt(2)^(t-1)-1))
            scrambles[i,j] = BTree(BTree(rbitsleft,seq.TInt(0)),BTree(rbitsright,seq.TInt(2)^(t-1)))
        end 
    end 
    RandomOwenScramble("Rand Owen Scramble: "*seq.name,seq,r,rngs,scrambles,t,t-seq.t,recipd)
end

RandomOwenScramble(seq::DigitalSeqB2G,r::Int64,rng::Xoshiro) = RandomOwenScramble(seq,r,reshape(spawn(rng,r*seq.s),r,seq.s))

RandomOwenScramble(seq::DigitalSeqB2G,r::Int64,seed::Int64) = RandomOwenScramble(seq,r,Xoshiro(seed))

RandomOwenScramble(seq::DigitalSeqB2G,r::Int64) = RandomOwenScramble(seq,r,Xoshiro())

RandomOwenScramble(seq::DigitalSeqB2G) = RandomOwenScramble(seq,1)

function GetScrambleScalar(xb::Union{UInt64,UInt128,BigInt},t::Int64,scramble::BTree,rng::Xoshiro,TInt::DataType)
    if scramble.xb === nothing # branch node, typeof(scramble.rbits) == Bool
        r1 = TInt(scramble.rbits)<<(t-1)
        b = Bool((xb>>(t-1))&1)
        onesmask = TInt(2)^(t-1)-1
        xbnext = xb&onesmask
        if ~b & (scramble.left === nothing)
            rbits = rand(rng,0:onesmask)
            scramble.left = BTree(rbits,xbnext)
            return r1+rbits
        elseif b & (scramble.right === nothing)
            rbits = rand(rng,0:onesmask)
            scramble.right = BTree(rbits,xbnext)
            return r1+rbits
        end
        scramble = ~b ? scramble.left : scramble.right
        return  r1 + GetScrambleScalar(xbnext,t-1,scramble,rng,TInt)
    elseif scramble.xb != xb # unseen leaf node
        ogsrbits,orsxb = scramble.rbits,scramble.xb
        b,ubit = nothing,nothing
        rmask = TInt(2)^t-1
        while true
            b,ubit,rbit = Bool((xb>>(t-1))&1),Bool((orsxb>>(t-1))&1),Bool((ogsrbits>>(t-1))&1)
            scramble.rbits,scramble.xb = rbit,nothing
            if ubit != b break end
            scramble = ~b ? scramble.left = BTree() : scramble.right = BTree()
            t -= 1
        end
        onesmask = TInt(2)^(t-1)-1
        newrbits = rand(rng,0:onesmask)
        scramble.left = ~b ? BTree(newrbits,xb&onesmask) : BTree(ogsrbits&onesmask,orsxb&onesmask)
        scramble.right = b ? BTree(newrbits,xb&onesmask) : BTree(ogsrbits&onesmask,orsxb&onesmask)
        rmask ⊻= onesmask
        return (ogsrbits&rmask) + newrbits
    else # scramble.xb == xb
        return scramble.rbits # seen leaf node 
    end
end

function OwenScramble(rds::RandomOwenScramble,xb::Union{Matrix{UInt64},Matrix{UInt128},Matrix{BigInt}},n::Int64)
    xrbs = [zeros(rds.seq.TInt,n,rds.seq.s) for j=1:rds.r]
    for i=1:n
        for l=1:rds.r
            for j=1:rds.seq.s
                xbij = xb[i,j]<<rds.tdiff
                b = Bool(xbij>>(rds.t-1)&1)
                scramble = ~b ? rds.scrambles[l,j].left : rds.scrambles[l,j].right
                xrbs[l][i,j] = xbij ⊻ GetScrambleScalar(xbij,rds.t,scramble,rds.rngs[l,j],rds.seq.TInt)
            end
        end 
    end
    xrbs
end 

NextRBinary(rds::RandomOwenScramble,n::Int64) = OwenScramble(rds,NextBinary(rds.seq,n),n)

FirstRLinearBinary(rds::RandomOwenScramble,m::Int64) = OwenScramble(rds,FirstLinearBinary(rds.seq,m),2^m)
