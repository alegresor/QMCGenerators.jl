include("qmc_generators.jl")
using .QMCGenerators

m = 5
C1 = [BigInt(2^i) for i=0:(m-1)]
C2 = [BigInt(1) for i=1:m]
for i in 2:m C2[i] = (C2[i-1] << 1) ⊻ C2[i-1] end
Cs = vcat(C1',C2')
seq = DigitalSeqB2G(2,Cs)
x = Next(seq,4)
display(x)

rseq = RandomDigitalShift(seq,3)
xr = NextR(rseq,4)
for j=1:size(xr,1)
    println("\nj=$j")
    display(xr[j])
end
println()

seq = DigitalSeqB2G(3,"nxmats/nx_b2_m30_s4_Cs.col")
rseq = RandomDigitalShift(seq,1)
xr = Next(rseq,8)
display(xr)
println()

seq = DigitalSeqB2G(2)
x = FirstLinear(seq,3)
t = Next(seq,8)
display(t); println(); display(x)
