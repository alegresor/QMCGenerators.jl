function plot_mc_vs_qmc(m = 16, r = 100, output = true)
    seed = 7
    n = 2^m
    s,mu = 7,-11.05684907978818
    colors = ["#9558B2","#389826","#4063D8"]
    rng = MersenneTwister(seed)
    xsets = [
        ["Pseudo-Random", [rand(rng,n,s) for k=1:r]],
        ["Digital Net (Quasi-Random)", NextR(RandomDigitalShift(DigitalSeqB2G(s),r,seed),n)],
        ["Lattice (Quasi-Random)", NextR(RandomShift(LatticeSeqB2(s),r,seed),n)]
    ]

    if ~output return end 

    f(x::Vector{Float64}) = π^(s/2)*cos(norm(quantile.(Normal(),x)/sqrt(2)));
    f(x::Matrix{Float64}) = map(i->f(x[i,:]),1:size(x,1))

    fig,ax = PyPlot.subplots(nrows=1,ncols=1,figsize=(10,6))
    for k=1:size(xsets,1)
        name,xs = xsets[k]
        ys = vcat(map(i->f(xs[i]),1:r)'...)
        muhats = cumsum(ys,dims=2); for i=1:r muhats[i,:] = muhats[i,:]./range(1,n) end 
        err = abs.(muhats.-mu)
        pows2 = 2 .^ (0:m)
        qlowerr = map(p2->quantile(err[:,p2],.35),pows2)
        qmid = map(p2->quantile(err[:,p2],.5),pows2)
        qhigherr = map(p2->quantile(err[:,p2],.65),pows2)
        ax.plot(pows2,qmid,color=colors[k],label=name)
        ax.fill_between(pows2,qhigherr,qlowerr,color=colors[k],alpha=.25)
    end
    ax.set_xlim([1,n]); ax.set_xlabel("samples"); ax.set_ylabel("error")
    ax.set_xscale("log",base=2); ax.set_yscale("log",base=10)
    ax.spines["top"].set_visible(false); ax.spines["right"].set_visible(false)
    ax.legend(loc="lower left",frameon=false)
    fig.savefig(joinpath(@__DIR__(),"../docs/src/assets/mc_vs_qmc.svg"),bbox_inches="tight",transparent=true)
    PyPlot.close()
end
