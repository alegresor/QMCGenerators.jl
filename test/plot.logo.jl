function plot_logo(output = true)
    ns = 2 .^ [2,4,6]
    colors = ["#9558B2","#389826","#4063D8"]
    @assert length(ns)==length(colors)
    n = maximum(ns)
    x = Next(RandomDigitalShift(DigitalSeqB2G(2),1,17),n)

    if ~output return end

    fig,ax = PyPlot.subplots(figsize=(8,8))
    global nmin = 1
    for i=1:length(ns)
        n_max = ns[i]
        xi = x[nmin:n_max,:]
        ax.scatter(xi[:,1],xi[:,2],color=colors[i],s=500,marker=".")
        global nmin = n_max+1
    end
    ax.spines["top"].set_visible(false); ax.spines["bottom"].set_visible(false); ax.spines["left"].set_visible(false); ax.spines["right"].set_visible(false) 
    ax.set_xlim([0,1]); ax.set_ylim([0,1]); ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect(1)
    for i=1:7 ax.axvline(x=i/8,color=colors[3],alpha=1); ax.axhline(y=i/8,color=colors[3],alpha=1) end 
    for i=1:3 ax.axvline(x=i/4,color=colors[2],alpha=1); ax.axhline(y=i/4,color=colors[2],alpha=1) end 
    for i=1:1 ax.axvline(x=i/2,color=colors[1],alpha=1); ax.axhline(y=i/2,color=colors[1],alpha=1) end 
    fig.savefig(joinpath(@__DIR__(),"../docs/src/assets/logo.svg"),bbox_inches="tight",transparent=true)
    PyPlot.close()
end 