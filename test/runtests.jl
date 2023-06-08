using QMCGenerators
using Test

@testset "LatticeSeqB2" begin 
    x = Next(LatticeSeqB2(3),8)
    @test x == vcat([
        [0.0,    0.0,    0.0],
        [0.5,    0.5,    0.5],
        [0.25,   0.75,   0.75],
        [0.75,   0.25,   0.25],
        [0.125,  0.375,  0.375],
        [0.625,  0.875,  0.875],
        [0.375,  0.125,  0.125],
        [0.875,  0.625,  0.625]]'...)

    x = FirstLinear(LatticeSeqB2(3),3)
    @test x == vcat([
        [0.0,    0.0,    0.0],
        [0.125,  0.375,  0.375],
        [0.25,   0.75,   0.75],
        [0.375,  0.125,  0.125],
        [0.5,    0.5,    0.5],
        [0.625,  0.875,  0.875],
        [0.75,   0.25,   0.25],
        [0.875,  0.625,  0.625]]'...)

    m,s,r = 3,3,2
    n = 2^m

    xrs = NextR(RandomShift(LatticeSeqB2(s),r),n)
    @test size(xrs,1)==r
    @test all(size(xrs[i])==(n,s) for i=1:r)

    xrs = FirstRLinear(RandomShift(LatticeSeqB2(s),r),m)
    @test size(xrs,1)==r
    @test all(size(xrs[i])==(n,s) for i=1:r)

    xr = Next(RandomShift(LatticeSeqB2(s)),n)
    @test size(xr)==(n,s)
    @test Next(RandomShift(LatticeSeqB2(s),1,7),n) == Next(RandomShift(LatticeSeqB2(s),1,7),n)

    rls = RandomShift(LatticeSeqB2(s,"exod2_base2_m20.txt",20))
    x = Next(rls,n)
    @test size(x)==(n,s)
end 

@testset "DigitalSeqB2G.jl" begin
    m = 5
    C1 = [BigInt(2^i) for i=0:(m-1)]
    C2 = [BigInt(1) for i=1:m]
    for i in 2:m C2[i] = (C2[i-1] << 1) ⊻ C2[i-1] end
    Cs = vcat(C1',C2')
    seq = DigitalSeqB2G(2,Cs)
    x = Next(seq,8)
    @test x == vcat([
        [0.0,    0.0],
        [0.5,    0.5],
        [0.75,   0.25],
        [0.25,   0.75],
        [0.375,  0.375],
        [0.875,  0.875],
        [0.625,  0.125],
        [0.125,  0.625]]'...)

    x = Next(DigitalSeqB2G(3),8)
    @test x == vcat([
        [0.0,    0.0,    0.0],
        [0.5,    0.5,    0.5],
        [0.75,   0.25,   0.25],
        [0.25,   0.75,   0.75],
        [0.375,  0.375,  0.625],
        [0.875,  0.875,  0.125],
        [0.625,  0.125,  0.875],
        [0.125,  0.625,  0.375]]'...)

    x = FirstLinear(DigitalSeqB2G(3),3)
    @test x == vcat([
        [0.0,    0.0,    0.0],
        [0.5,    0.5,    0.5],
        [0.25,   0.75,   0.75],
        [0.75,   0.25,   0.25],
        [0.125,  0.625,  0.375],
        [0.625,  0.125,  0.875],
        [0.375,  0.375,  0.625],
        [0.875,  0.875,  0.125]]'...)

    m,s,r = 3,3,2
    n = 2^m

    xrs = NextR(RandomDigitalShift(DigitalSeqB2G(s),r),n)
    @test size(xrs,1)==r
    @test all(size(xrs[i])==(n,s) for i=1:r)

    xrs = FirstRLinear(RandomDigitalShift(DigitalSeqB2G(s),r),m)
    @test size(xrs,1)==r
    @test all(size(xrs[i])==(n,s) for i=1:r)

    xr = Next(RandomDigitalShift(DigitalSeqB2G(s)),n)
    @test size(xr)==(n,s)
    @test Next(RandomDigitalShift(DigitalSeqB2G(s),1,7),n) == Next(RandomDigitalShift(DigitalSeqB2G(s),1,7),n)

    rls = RandomDigitalShift(DigitalSeqB2G(3,"sobol_mini.col"))
    rls = RandomDigitalShift(DigitalSeqB2G(s,"nxmats/nx_b2_m30_s10_Cs.col"))
    rls = RandomDigitalShift(DigitalSeqB2G(s,"sobolmats/sobol_alpha3_Bs64.col"))
    x = Next(rls,n)
    @test size(x)==(n,s)
end
