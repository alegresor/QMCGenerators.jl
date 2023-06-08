include("../../../src/QMCGenerators.jl")
using .QMCGenerators
using PyPlot

ns = 2 .^ [2,4,6]
colors = ["#9558B2","#389826","#4063D8"]
@assert length(ns)==length(colors)
n = maximum(ns)
x = Next(RandomDigitalShift(DigitalSeqB2G(2),1,17),n)

fig,ax = PyPlot.subplots(figsize=(8,8))
nmin = 1
for i=1:length(ns)
    n_max = ns[i]
    xi = x[nmin:n_max,:]
    ax.scatter(xi[:,1],xi[:,2],color=colors[i],s=500,marker=".")
    global nmin = n_max+1
end
ax.spines["top"].set_visible(false); ax.spines["bottom"].set_visible(false); ax.spines["left"].set_visible(false); ax.spines["right"].set_visible(false) 
ax.set_xlim([0,1]); ax.set_ylim([0,1]); ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect(1)
for i=0:8 ax.axvline(x=i/8,color="k",alpha=.25); ax.axhline(y=i/8,color="k",alpha=.25) end 
for i=0:4 ax.axvline(x=i/4,color="k",alpha=.5); ax.axhline(y=i/4,color="k",alpha=.5) end 
for i=0:2 ax.axvline(x=i/2,color="k",alpha=.75); ax.axhline(y=i/2,color="k",alpha=.75) end 
fig.savefig(joinpath(@__DIR__(),"logo.svg"),bbox_inches="tight",transparent=true)
