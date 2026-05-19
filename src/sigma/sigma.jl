module Sigma
using Cuba
using JLD2, CSV

using Printf, LinearAlgebra, DataFrames
using ..CompositeGrids
using ..ElectronGas
using ..MCIntegration
using ..Lehmann

using ..FeynmanDiagram
import ..FeynmanDiagram.FrontEnds: Filter, NoHartree, NoFock, DirectOnly, Wirreducible, Girreducible, NoBubble, Proper
import ..FeynmanDiagram.Parquet: DiagPara, SigmaDiag
using ..Measurements

using ..UEG
using ..Propagator
import ..Propagator: LeafStateAD
using ..Diagram

@inline function phase(varT, extT, l, β)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return exp(1im * π * (2l + 1) / β * (tout - tin))
end

include("parquetAD.jl")
include("sigmaKW.jl")
include("sigma_dk.jl")
include("sigmaKW_Reweight.jl")
include("sigmadk_Reweight.jl")
# include("sigmaCuba.jl")
# include("sigmaVegas.jl")

function MC_Clib(para; kgrid=[para.kF,], ngrid=[0], neval=1e6, reweight_goal=nothing,
    # spinPolarPara::Float64=0.0, # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
    filename::Union{String,Nothing}=nothing, partition=UEG.partition(para.order),
    isLayered2D=false, # whether to use the screened Coulomb interaction in 2D or not 
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"), verbose=-1
)
    @assert para.spin == 2 "Only spin-unpolarized case is supported for compiled C library"
    kF = para.kF
    neighbor = UEG.neighbor(partition)

    if isLayered2D
        @assert (para.dim == 2) "Only 2D systems supports the tanh screened Coulomb interaction"
    end

    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            reweight_factor = 2.0^(2order + sOrder + vOrder - 2)
            if (order, sOrder, vOrder) == (1, 0, 0)
                reweight_factor = 4.0
            end
            push!(reweight_goal, reweight_factor)
        end
        push!(reweight_goal, 4.0)
    end

    diaginfo = Sigma.diagram_loadinfo(para, partition, root_dir=root_dir)
    sigma, result = Sigma.ParquetAD_Clib(para, diaginfo;
        root_dir=root_dir, isLayered2D=isLayered2D,
        neighbor=neighbor, reweight_goal=reweight_goal,
        kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread, print=verbose)

    if isnothing(sigma) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (ngrid, kgrid, sigma)
            end
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            @printf("%10s  %10s   %10s   %10s   %10s \n", "q/kF", "real(avg)", "err", "imag(avg)", "err")
            r, i = real(sigma[key]), imag(sigma[key])
            for (in, n) in enumerate(ngrid)
                println("n = $n")
                for (iq, q) in enumerate(kgrid)
                    @printf("%10.6f  %10.6f ± %10.6f   %10.6f ± %10.6f\n", q[1] / kF, r[in, iq].val, r[in, iq].err, i[in, iq].val, i[in, iq].err)
                end
            end
        end
    end
    return sigma, result
end

function MC(para; kgrid=[para.kF,], ngrid=[0], neval=1e6, reweight_goal=nothing,
    # spinPolarPara::Float64=0.0, # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
    filename::Union{String,Nothing}=nothing, partition=UEG.partition(para.order),
    isLayered2D=false, # whether to use the screened Coulomb interaction in 2D or not 
    filter=[NoHartree], extK=nothing, optimize_level=1, verbose=-1
)
    kF = para.kF
    neighbor = UEG.neighbor(partition)

    if isLayered2D
        @assert (para.dim == 2) "Only 2D systems supports the tanh screened Coulomb interaction"
    end

    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            reweight_factor = 2.0^(2order + sOrder + vOrder - 2)
            if (order, sOrder, vOrder) == (1, 0, 0)
                reweight_factor = 4.0
            end
            push!(reweight_goal, reweight_factor)
        end
        push!(reweight_goal, 4.0)
    end

    diagram = Diagram.diagram_parquet_noresponse(:sigma, para, partition, filter=filter, extK=extK, optimize_level=optimize_level)
    sigma, result = Sigma.ParquetAD(para, diagram;
        isLayered2D=isLayered2D, print=verbose,
        neighbor=neighbor, reweight_goal=reweight_goal,
        kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread)

    if isnothing(sigma) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (ngrid, kgrid, sigma)
            end
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            @printf("%10s  %10s   %10s   %10s   %10s \n", "q/kF", "real(avg)", "err", "imag(avg)", "err")
            r, i = real(sigma[key]), imag(sigma[key])
            for (in, n) in enumerate(ngrid)
                println("n = $n")
                for (iq, q) in enumerate(kgrid)
                    @printf("%10.6f  %10.6f ± %10.6f   %10.6f ± %10.6f\n", q[1] / kF, r[in, iq].val, r[in, iq].err, i[in, iq].val, i[in, iq].err)
                end
            end
        end
    end
    return sigma, result
end

function diagram_loadinfo(paramc::ParaMC, _partition::Vector{T};
    filter=[NoHartree], transferLoop=nothing,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"), filename="extvars_sigma.jld2"
) where {T}
    diagpara = Vector{DiagPara}()
    extT_labels = Vector{Vector{Int}}[]

    fname = joinpath(root_dir, filename)
    jldopen(fname, "r") do f
        for p in _partition
            key_str = join(string.(p))
            if key_str in keys(f)
                extT = f[key_str][1]
                push!(diagpara, Diagram.diagPara(SigmaDiag, paramc.isDynamic, p[1], paramc.spin, filter, transferLoop))
                push!(extT_labels, extT)
            else
                error("$(key_str) not found in $(fname)")
            end
        end
    end
    return (_partition, diagpara, extT_labels)
end

# include("source_codeParquetAD/Cwrapper_sigma_ParquetAD.jl")
# include("source_codeParquetAD/Cwrapper_sigmadk_ParquetAD.jl")
include("source_codeGV/Cwrapper_sigma_GV.jl")
include("source_codeGV/Cwrapper_sigmadk_GV.jl")

const evalfuncParquetAD_sigma_map = Dict{Tuple{Int,Int,Int},Function}()
const evalfuncParquetAD_sigmadk_map = Dict{Tuple{Int,Int,Int,Int},Function}()

for sym in names(@__MODULE__; all=true)
    s_str = String(sym)

    m_sigma = match(r"^eval_sigma_ParquetAD(\d)(\d)(\d)!$", s_str)
    if m_sigma !== nothing
        o, i, j = parse.(Int, m_sigma.captures)
        evalfuncParquetAD_sigma_map[(o, i, j)] = getfield(@__MODULE__, sym)
        continue
    end

    m_sigmadk = match(r"^eval_sigmadk_ParquetAD(\d)(\d)(\d)(\d)!$", s_str)
    if m_sigmadk !== nothing
        o, i, j, k = parse.(Int, m_sigmadk.captures)
        evalfuncParquetAD_sigmadk_map[(o, i, j, k)] = getfield(@__MODULE__, sym)
    end
end

end
