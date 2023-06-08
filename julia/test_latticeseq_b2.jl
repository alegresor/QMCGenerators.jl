include("qmc_generators.jl")
using .QMCGenerators

ls = LatticeSeqB2(3)
x = Next(ls,4)
display(x)

rls = RandomShift(LatticeSeqB2(3),2,5)
xr = NextR(rls,4)
for j=1:size(xr,1)
    println("\nj=$j")
    display(xr[j])
end 
println()

rls = RandomShift(LatticeSeqB2(3,"exod2_base2_m20.txt",20))
x = Next(rls,8)
display(x)
println()

rls = RandomShift(LatticeSeqB2(3),2)
xr = FirstRLinear(rls,2)
for j=1:size(xr,1)
    println("\nj=$j")
    display(xr[j])
end 
println()

rls = RandomShift(LatticeSeqB2(3,"exod2_base2_m20.txt",20))
x = FirstLinear(rls,2)
display(x)