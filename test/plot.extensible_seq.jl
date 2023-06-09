function plot_extensible_seq(ns = 2 .^ [6,7,8], output = true)
    @assert length(ns) == 3
    seed = 7
    n = maximum(ns)
    colors = ["#9558B2","#389826","#4063D8"]
    xsets = [
        ["Pseudo-Random", rand(MersenneTwister(seed),n,2)],
        ["Digital Net (Quasi-Random)", Next(RandomDigitalShift(DigitalSeqB2G(2),1,seed),n)],
        ["Lattice (Quasi-Random)", Next(RandomShift(LatticeSeqB2(2),1,seed),n)]
    ]

    if ~output return end 

    ncols = size(xsets,1)
    fig,ax = PyPlot.subplots(nrows=1,ncols=ncols,figsize=(4*ncols,4))
    nmin = NaN
    for i=1:ncols
        name,xfull = xsets[i]
        local nmin = 1
        for j=1:size(ns,1)
            nmax = ns[j]
            x = xfull[nmin:nmax,:]
            ax[i].scatter(x[:,1],x[:,2],s=10,color=colors[j])
            nmin = nmax+1
        end
        ax[i].set_title(name)
        ax[i].set_xlim([0,1]); ax[i].set_ylim([0,1])
        ax[i].set_xticks([0,1]); ax[i].set_yticks([0,1])
        ax[i].set_aspect(1)
    end 
    fig.savefig(joinpath(@__DIR__(),"../docs/src/assets/extensible_seq.svg"),bbox_inches="tight",transparent=false) 
    PyPlot.close()
end 