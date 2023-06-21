JULIA4LOGOCOLORS = ["#4063D8","#389826","#9558B2","#CB3C33"]

function qmcscatter!(ax::CairoMakie.Axis,x::Matrix{Float64},nvec::Vector{Int64})
    lnvec = length(nvec)
    @assert size(x,1) >= maximum(nvec) && size(x,2) == 2
    @assert lnvec <= (length(JULIA4LOGOCOLORS)+1)
    for i=1:lnvec-1
        nmin,nmax = i == 1 ? nvec[1] : nvec[i]+1,nvec[i+1]
        CairoMakie.scatter!(ax,x[nmin:nmax,1],x[nmin:nmax,2],color=JULIA4LOGOCOLORS[i],label="$nmin:$nmax")
    end
end

function qmcscatter!(xs::Vector{Matrix{Float64}},nvec::Vector{Int64},dvec::Matrix{Int64}=[1 2;])
    @assert all(size(xs[j],2) >= maximum(dvec) for j=1:length(xs))
    nrows,ncols = size(dvec,1),length(xs)
    fig = CairoMakie.Figure(backgroundcolor=:white, resolution=(400*ncols,400*nrows))
    ax = NaN
    for i=1:nrows
        for j=1:ncols
            dhoriz,dvert = dvec[i,1],dvec[i,2]
            x = xs[j][:,[dhoriz,dvert]]
            ax = CairoMakie.Axis(fig[i,j],
                xlabel = latexstring("\$X_{i\\;$dhoriz}\$"),
                ylabel = latexstring( "\$X_{i\\;$dvert}\$"),
                xticks = ([0,1/4,1/2,3/4,1],[L"$0$",L"$1/4$",L"$1/2$",L"$3/4$",L"$1$"]),
                yticks = ([0,1/4, 1/2,3/4,1],[L"$0$",L"$1/4$",L"$1/2$",L"$3/4$",L"$1$"]),
                aspect = 1)
            if i==1 && ncols>1 ax.title = "Randomization $j" end 
            CairoMakie.limits!(ax,[0,1],[0,1])
            qmcscatter!(ax,x,nvec)
        end
    end
    if length(nvec)>2 fig[nrows+1,:] = CairoMakie.Legend(fig,ax,"Index",framevisible=false,orientation=:horizontal) end
    return fig
end

qmcscatter!(xs::Vector{Matrix{Float64}},nvec::Int64,args...) = qmcscatter!(xs,[1,nvec],args...)

function qmcscatter!(seq::Union{DigitalSeqB2G,LatticeSeqB2},nvec::Union{Int64,Vector{Int64}},args...)
    @assert seq.k == -1
    xs = [Next(seq,maximum(nvec))]
    Reset!(seq)
    qmcscatter!(xs,nvec,args...)
end

function qmcscatter!(rseq::Union{RandomDigitalShift,RandomShift},nvec::Union{Int64,Vector{Int64}},args...)
    if typeof(rseq) == RandomDigitalShift @assert rseq.ds.k == -1 end 
    if typeof(rseq) == RandomShift @assert rseq.ls.k == -1 end 
    xs = NextR(rseq,maximum(nvec))
    Reset!(rseq)
    qmcscatter!(xs,nvec,args...)
end