mutable struct LinearMatrixScramble
    s::Int64 
    Csrlms::Matrix{BigInt}
    m::Int64
    t::Int64 
end

function LinearMatrixScramble(s::Int64,Cs::Matrix{BigInt},rng::Xoshiro) # Matousek
    @assert s <= size(Cs,1)
    m = size(Cs,2)
    tog = maximum(ndigits.(Cs[1,:],base=2))
    t = max(53,tog)
    Csr = bitreverse.(Cs[1:s,:],tog)
    Csrlms = zeros(BigInt,s,m)
    for j=1:s
        for k=0:t-1
            k1 = min(k,tog)
            delem = BigInt(1)<<(tog-k1-1)
            ndelem = rand(rng,0:(BigInt(1)<<k1)-1)<<(tog-k1)
            u = delem+ndelem
            for l=1:m
                v = u&Csr[j,l]
                setbits = 0 
                for i=0:tog setbits += (v>>i)&1 end 
                b = setbits%2
                Csrlms[j,l] += BigInt(b)<<(t-k-1)
            end
        end
    end
    LinearMatrixScramble(s,Csrlms,m,t)
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
    name::String
    s::Int64 # dimension
    Csr::Matrix{BigInt} # matrix of size at least (s,m)
    m::Int64 # number of columns, can generate 2^m points
    t::Int64 # maximum number of bits in an element of Cs
    alpha::Float64 # t/m, the order of the net 
    n::Int64 # maximum number of supported points
    recipd::Union{BigFloat,Float64} # multiplication factor
    k::Int64 # index in the sequence
    cur::Array{BigInt}
end

function DigitalSeqB2G(s::Int64,Cs::Matrix{BigInt})
    @assert s <= size(Cs,1)
    m = size(Cs,2)
    t = maximum(ndigits.(Cs[1,:],base=2))
    Csr = bitreverse.(Cs[1:s,:],t)
    alpha = t/m 
    n = 2^m
    recipd = t>53 ? BigFloat(2)^(-t) : Float64(2)^(-t)
    cur = zeros(BigInt,s)
    DigitalSeqB2G("Digital Seq B2",s,Csr,m,t,alpha,n,recipd,-1,cur)
end

DigitalSeqB2G(s::Int64,path::String) = DigitalSeqB2G(s,readdlm(download(joinpath("https://bitbucket.org/dnuyens/qmc-generators/raw/cb0f2fb10fa9c9f2665e41419097781b611daa1e/DIGSEQ/",path)),BigInt))

DigitalSeqB2G(s::Int64) = DigitalSeqB2G(s,DEFAULT_DIGITALSEQB2G_GMATRIX)

function DigitalSeqB2G(rlms::LinearMatrixScramble)
    alpha = rlms.t/rlms.m
    n = 2^rlms.m
    recipd = rlms.t>53 ? BigFloat(2)^(-rlms.t) : Float64(2)^(-rlms.t)
    cur = zeros(BigInt,rlms.s)
    DigitalSeqB2G("LMS Digital Seq B2",rlms.s,rlms.Csrlms,rlms.m,rlms.t,alpha,n,recipd,-1,cur)
end 

function Reset!(seq::DigitalSeqB2G) 
    seq.k = -1 
    seq.cur = zeros(BigInt,seq.s)
    return 
end 

function NextBinary(seq::DigitalSeqB2G,n::Int64)
    xb = zeros(BigInt,n,seq.s)
    for i=1:n
        seq.k += 1
        seq.k == seq.n && throw(DomainError(seq.k,"already generated maximum number of points"))
        if seq.k ==0 
            xb[i,:] = seq.cur # zeros
            continue 
        end 
        ctz = ndigits(((seq.k ⊻ (seq.k-1))+1) >> 1, base=2) - 1
        seq.cur .⊻= seq.Csr[:,(ctz+1)]
        xb[i,:] = seq.cur
    end
    xb
end

Next(seq::DigitalSeqB2G,n::Int64) = BinaryToFloat64(NextBinary(seq,n),seq)

function FirstLinearBinary(seq::DigitalSeqB2G,m::Int64)
    @assert seq.k==-1
    n = 2^m
    xb = NextBinary(seq,n)
    gcs = map(i->i⊻(i>>1),0:n-1)
    xbl = zeros(BigInt,size(xb))
    for i=1:n xbl[gcs[i]+1,:] .= xb[i,:] end 
    Reset!(seq)
    xbl
end

FirstLinear(seq::DigitalSeqB2G,m::Int64) = BinaryToFloat64(FirstLinearBinary(seq,m),seq)

mutable struct RandomDigitalShift
    name::String
    seq::DigitalSeqB2G
    r::Int64
    rshifts::Matrix{BigInt}
    t::Int64 # number of bits in shifted integers
    tdiff::Int64
    recipd::Union{BigFloat,Float64} # multiplication factor
end 

function RandomDigitalShift(seq::DigitalSeqB2G,r::Int64,rng::Xoshiro)
    t = max(53,seq.t)
    recipd = t>53 ? BigFloat(2)^(-t) : Float64(2)^(-t)
    rshifts = rand(rng,0:(BigInt(2)^t-1),r,seq.s)
    RandomDigitalShift("Rand DShift: "*seq.name,seq,r,rshifts,t,t-seq.t,recipd)
end

RandomDigitalShift(seq::DigitalSeqB2G,r::Int64,seed::Int64) = RandomDigitalShift(seq,r,Xoshiro(seed))

RandomDigitalShift(seq::DigitalSeqB2G,r::Int64) = RandomDigitalShift(seq,r,Xoshiro())

RandomDigitalShift(seq::DigitalSeqB2G) = RandomDigitalShift(seq,1)

DigitalShifts(xb::Matrix{BigInt},rds::RandomDigitalShift) = [rds.rshifts[i,:]' .⊻ xb for i=1:rds.r]

NextRBinary(rds::RandomDigitalShift,n::Int64) = DigitalShifts(NextBinary(rds.seq,n).<<rds.tdiff,rds)

FirstRLinearBinary(rds::RandomDigitalShift,m::Int64) = DigitalShifts(FirstLinearBinary(rds.seq,m).<<rds.tdiff,rds)

mutable struct BTree
    rbits::Union{Nothing,BigInt,Bool}
    xb::Union{Nothing,BigInt}
    left::Union{Nothing,BTree}
    right::Union{Nothing,BTree}
end
BTree(left::BTree,right::BTree) = BTree(nothing,nothing,left,right)
BTree(rbits::BigInt,xb::BigInt) = BTree(rbits,xb,nothing,nothing)
BTree() = BTree(nothing,nothing,nothing,nothing)

mutable struct RandomOwenScramble
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
            r1 = BigInt(rand(rngs[i,j],Bool))<<(t-1)
            rbitsleft,rbitsright = r1 + rand(rngs[i,j],0:(BigInt(2)^(t-1)-1)),r1 + rand(rngs[i,j],0:(BigInt(2)^(t-1)-1))
            scrambles[i,j] = BTree(BTree(rbitsleft,BigInt(0)),BTree(rbitsright,BigInt(2)^(t-1)))
        end 
    end 
    RandomOwenScramble("Rand Owen Scramble: "*seq.name,seq,r,rngs,scrambles,t,t-seq.t,recipd)
end

RandomOwenScramble(seq::DigitalSeqB2G,r::Int64,rng::Xoshiro) = RandomOwenScramble(seq,r,reshape(spawn(rng,r*seq.s),r,seq.s))

RandomOwenScramble(seq::DigitalSeqB2G,r::Int64,seed::Int64) = RandomOwenScramble(seq,r,Xoshiro(seed))

RandomOwenScramble(seq::DigitalSeqB2G,r::Int64) = RandomOwenScramble(seq,r,Xoshiro())

RandomOwenScramble(seq::DigitalSeqB2G) = RandomOwenScramble(seq,1)

function GetScrambleScalar(xb::BigInt,t::Int64,scramble::BTree,rng::Xoshiro)
    if scramble.xb === nothing # branch node, typeof(scramble.rbits) == Bool
        r1 = BigInt(scramble.rbits)<<(t-1)
        b = Bool((xb>>(t-1))&1)
        onesmask = BigInt(2)^(t-1)-1
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
        return  r1 + GetScrambleScalar(xbnext,t-1,scramble,rng)
    elseif scramble.xb != xb # unseen leaf node
        ogsrbits,orsxb = scramble.rbits,scramble.xb
        b,ubit = nothing,nothing
        rmask = BigInt(2)^t-1
        while true
            b,ubit,rbit = Bool((xb>>(t-1))&1),Bool((orsxb>>(t-1))&1),Bool((ogsrbits>>(t-1))&1)
            scramble.rbits,scramble.xb = rbit,nothing
            if ubit != b break end
            scramble = ~b ? scramble.left = BTree() : scramble.right = BTree()
            t -= 1
        end
        onesmask = BigInt(2)^(t-1)-1
        newrbits = rand(rng,0:onesmask)
        scramble.left = ~b ? BTree(newrbits,xb&onesmask) : BTree(ogsrbits&onesmask,orsxb&onesmask)
        scramble.right = b ? BTree(newrbits,xb&onesmask) : BTree(ogsrbits&onesmask,orsxb&onesmask)
        rmask ⊻= onesmask
        return (ogsrbits&rmask) + newrbits
    else # scramble.xb == xb
        return scramble.rbits # seen leaf node 
    end
end

function OwenScramble(rds::RandomOwenScramble,xb::Matrix{BigInt},n::Int64)
    xrbs = [zeros(BigInt,n,rds.seq.s) for j=1:rds.r]
    for i=1:n
        for l=1:rds.r
            for j=1:rds.seq.s
                xbij = xb[i,j]<<rds.tdiff
                b = Bool(xbij>>(rds.t-1)&1)
                scramble = ~b ? rds.scrambles[l,j].left : rds.scrambles[l,j].right
                xrbs[l][i,j] = xbij ⊻ GetScrambleScalar(xbij,rds.t,scramble,rds.rngs[l,j])
            end
        end 
    end
    xrbs
end 

NextRBinary(rds::RandomOwenScramble,n::Int64) = OwenScramble(rds,NextBinary(rds.seq,n),n)

FirstRLinearBinary(rds::RandomOwenScramble,m::Int64) = OwenScramble(rds,FirstLinearBinary(rds.seq,m),2^m)
