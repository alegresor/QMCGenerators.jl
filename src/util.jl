function qmcscatter!(ax::CairoMakie.Axis,x::Matrix{Float64},nvec::Vector{Int64})
    @assert size(x,1)>=maximum(nvec) && size(x,2)==2
    nprev = 1
    for i in eachindex(nvec)
        n = nvec[i]
        CairoMakie.scatter!(ax,x[nprev:n,1],x[nprev:n,2])
        nprev = n+1
    end
end

function qmcscatter!(xs::Vector{Matrix{Float64}},nvec::Vector{Int64},dvec::Matrix{Int64}=[1 2;])
    @assert all(size(xs[j],2) >= maximum(dvec) for j=1:length(xs))
    nrows,ncols = size(dvec,1),length(xs)
    fig = CairoMakie.Figure(backgroundcolor=:white, resolution=(400*ncols,400*nrows))
    for i=1:nrows
        for j=1:ncols
            dhoriz,dvert = dvec[i,1],dvec[i,2]
            x = xs[j][:,[dhoriz,dvert]]
            ax = CairoMakie.Axis(fig[i,j],
                xlabel = latexstring("\$X_{i$dhoriz}\$"),
                ylabel = latexstring( "\$X_{i$dvert}\$"),
                xticks = ([0,1/4,1/2,3/4,1],[L"$0$",L"$1/4$",L"$1/2$",L"$3/4$",L"$1$"]),
                yticks = ([0,1/4, 1/2,3/4,1],[L"$0$",L"$1/4$",L"$1/2$",L"$3/4$",L"$1$"]),
                aspect = 1)
            CairoMakie.limits!(ax,[0,1],[0,1])
            qmcscatter!(ax,x,nvec)
        end 
    end
    return fig
end

qmcscatter!(xs::Vector{Matrix{Float64}},nvec::Int64,args...) = qmcscatter!(xs,[nvec],args...)

function qmcscatter!(seq::Union{DigitalSeqB2G,LatticeSeqB2},nvec::Union{Int64,Vector{Int64}},args...)
    Reset!(seq)
    xs = [Next(seq,maximum(nvec))]
    Reset!(seq)
    qmcscatter!(xs,nvec,args...)
end

function qmcscatter!(rseq::Union{RandomDigitalShift,RandomShift},nvec::Union{Int64,Vector{Int64}},args...)
    Reset!(rseq)
    xs = NextR(rseq,maximum(nvec))
    Reset!(rseq)
    qmcscatter!(xs,nvec,args...)
end