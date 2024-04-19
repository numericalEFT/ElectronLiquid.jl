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
    diag_generator::Symbol=:parquet,
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

    if diag_generator == :Parquet
        diagram = Diagram.diagram_parquet_noresponse(:sigma, para, partition, filter=filter, extK=extK, optimize_level=optimize_level)
    elseif diag_generator == :GV
        diagram = Diagram.diagram_GV_noresponse(:sigma, para, partition, filter=filter, optimize_level=optimize_level)
    else
        error("Unknown diag_generator: $diag_generator")
    end
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

include("source_codeParquetAD/Cwrapper_sigma_ParquetAD.jl")
include("source_codeParquetAD/Cwrapper_sigmadk_ParquetAD.jl")
# include("source_codeGV/Cwrapper_sigmaGV.jl")
# include("source_codeGV/Cwrapper_sigmaParquetAD_dk.jl")

const evalfuncParquetAD_sigma_map = Dict(
    (1, 0, 0) => eval_sigma_ParquetAD100!,
    (1, 0, 1) => eval_sigma_ParquetAD101!,
    (1, 0, 2) => eval_sigma_ParquetAD102!,
    (1, 0, 3) => eval_sigma_ParquetAD103!,
    (1, 0, 4) => eval_sigma_ParquetAD104!,
    (1, 0, 5) => eval_sigma_ParquetAD105!,
    (1, 1, 0) => eval_sigma_ParquetAD110!,
    (1, 1, 1) => eval_sigma_ParquetAD111!,
    (1, 1, 2) => eval_sigma_ParquetAD112!,
    (1, 1, 3) => eval_sigma_ParquetAD113!,
    (1, 1, 4) => eval_sigma_ParquetAD114!,
    (1, 2, 0) => eval_sigma_ParquetAD120!,
    (1, 2, 1) => eval_sigma_ParquetAD121!,
    (1, 2, 2) => eval_sigma_ParquetAD122!,
    (1, 2, 3) => eval_sigma_ParquetAD123!,
    (1, 3, 0) => eval_sigma_ParquetAD130!,
    (1, 3, 1) => eval_sigma_ParquetAD131!,
    (1, 3, 2) => eval_sigma_ParquetAD132!,
    (1, 4, 0) => eval_sigma_ParquetAD140!,
    (1, 4, 1) => eval_sigma_ParquetAD141!,
    (1, 5, 0) => eval_sigma_ParquetAD150!,
    (2, 0, 0) => eval_sigma_ParquetAD200!,
    (2, 0, 1) => eval_sigma_ParquetAD201!,
    (2, 0, 2) => eval_sigma_ParquetAD202!,
    (2, 0, 3) => eval_sigma_ParquetAD203!,
    (2, 0, 4) => eval_sigma_ParquetAD204!,
    (2, 1, 0) => eval_sigma_ParquetAD210!,
    (2, 1, 1) => eval_sigma_ParquetAD211!,
    (2, 1, 2) => eval_sigma_ParquetAD212!,
    (2, 1, 3) => eval_sigma_ParquetAD213!,
    (2, 2, 0) => eval_sigma_ParquetAD220!,
    (2, 2, 1) => eval_sigma_ParquetAD221!,
    (2, 2, 2) => eval_sigma_ParquetAD222!,
    (2, 3, 0) => eval_sigma_ParquetAD230!,
    (2, 3, 1) => eval_sigma_ParquetAD231!,
    (2, 4, 0) => eval_sigma_ParquetAD240!,
    (3, 0, 0) => eval_sigma_ParquetAD300!,
    (3, 0, 1) => eval_sigma_ParquetAD301!,
    (3, 0, 2) => eval_sigma_ParquetAD302!,
    (3, 0, 3) => eval_sigma_ParquetAD303!,
    (3, 1, 0) => eval_sigma_ParquetAD310!,
    (3, 1, 1) => eval_sigma_ParquetAD311!,
    (3, 1, 2) => eval_sigma_ParquetAD312!,
    (3, 2, 0) => eval_sigma_ParquetAD320!,
    (3, 2, 1) => eval_sigma_ParquetAD321!,
    (3, 3, 0) => eval_sigma_ParquetAD330!,
    (4, 0, 0) => eval_sigma_ParquetAD400!,
    (4, 0, 1) => eval_sigma_ParquetAD401!,
    (4, 0, 2) => eval_sigma_ParquetAD402!,
    (4, 1, 0) => eval_sigma_ParquetAD410!,
    (4, 1, 1) => eval_sigma_ParquetAD411!,
    (4, 2, 0) => eval_sigma_ParquetAD420!,
    (5, 0, 0) => eval_sigma_ParquetAD500!,
    (5, 0, 1) => eval_sigma_ParquetAD501!,
    (5, 1, 0) => eval_sigma_ParquetAD510!,
    (6, 0, 0) => eval_sigma_ParquetAD600!
)
const evalfuncParquetAD_sigmadk_map = Dict(
    (1, 0, 0, 1) => eval_sigmadk_ParquetAD1001!,
    (1, 0, 1, 1) => eval_sigmadk_ParquetAD1011!,
    (1, 0, 2, 1) => eval_sigmadk_ParquetAD1021!,
    (1, 0, 3, 1) => eval_sigmadk_ParquetAD1031!,
    (1, 0, 4, 1) => eval_sigmadk_ParquetAD1041!,
    (1, 0, 5, 1) => eval_sigmadk_ParquetAD1051!,
    (1, 1, 0, 1) => eval_sigmadk_ParquetAD1101!,
    (1, 1, 1, 1) => eval_sigmadk_ParquetAD1111!,
    (1, 1, 2, 1) => eval_sigmadk_ParquetAD1121!,
    (1, 1, 3, 1) => eval_sigmadk_ParquetAD1131!,
    (1, 1, 4, 1) => eval_sigmadk_ParquetAD1141!,
    (1, 2, 0, 1) => eval_sigmadk_ParquetAD1201!,
    (1, 2, 1, 1) => eval_sigmadk_ParquetAD1211!,
    (1, 2, 2, 1) => eval_sigmadk_ParquetAD1221!,
    (1, 2, 3, 1) => eval_sigmadk_ParquetAD1231!,
    (1, 3, 0, 1) => eval_sigmadk_ParquetAD1301!,
    (1, 3, 1, 1) => eval_sigmadk_ParquetAD1311!,
    (1, 3, 2, 1) => eval_sigmadk_ParquetAD1321!,
    (1, 4, 0, 1) => eval_sigmadk_ParquetAD1401!,
    (1, 4, 1, 1) => eval_sigmadk_ParquetAD1411!,
    (1, 5, 0, 1) => eval_sigmadk_ParquetAD1501!,
    (2, 0, 0, 1) => eval_sigmadk_ParquetAD2001!,
    (2, 0, 1, 1) => eval_sigmadk_ParquetAD2011!,
    (2, 0, 2, 1) => eval_sigmadk_ParquetAD2021!,
    (2, 0, 3, 1) => eval_sigmadk_ParquetAD2031!,
    (2, 0, 4, 1) => eval_sigmadk_ParquetAD2041!,
    (2, 1, 0, 1) => eval_sigmadk_ParquetAD2101!,
    (2, 1, 1, 1) => eval_sigmadk_ParquetAD2111!,
    (2, 1, 2, 1) => eval_sigmadk_ParquetAD2121!,
    (2, 1, 3, 1) => eval_sigmadk_ParquetAD2131!,
    (2, 2, 0, 1) => eval_sigmadk_ParquetAD2201!,
    (2, 2, 1, 1) => eval_sigmadk_ParquetAD2211!,
    (2, 2, 2, 1) => eval_sigmadk_ParquetAD2221!,
    (2, 3, 0, 1) => eval_sigmadk_ParquetAD2301!,
    (2, 3, 1, 1) => eval_sigmadk_ParquetAD2311!,
    (2, 4, 0, 1) => eval_sigmadk_ParquetAD2401!,
    (3, 0, 0, 1) => eval_sigmadk_ParquetAD3001!,
    (3, 0, 1, 1) => eval_sigmadk_ParquetAD3011!,
    (3, 0, 2, 1) => eval_sigmadk_ParquetAD3021!,
    (3, 0, 3, 1) => eval_sigmadk_ParquetAD3031!,
    (3, 1, 0, 1) => eval_sigmadk_ParquetAD3101!,
    (3, 1, 1, 1) => eval_sigmadk_ParquetAD3111!,
    (3, 1, 2, 1) => eval_sigmadk_ParquetAD3121!,
    (3, 2, 0, 1) => eval_sigmadk_ParquetAD3201!,
    (3, 2, 1, 1) => eval_sigmadk_ParquetAD3211!,
    (3, 3, 0, 1) => eval_sigmadk_ParquetAD3301!,
    (4, 0, 0, 1) => eval_sigmadk_ParquetAD4001!,
    (4, 0, 1, 1) => eval_sigmadk_ParquetAD4011!,
    (4, 0, 2, 1) => eval_sigmadk_ParquetAD4021!,
    (4, 1, 0, 1) => eval_sigmadk_ParquetAD4101!,
    (4, 1, 1, 1) => eval_sigmadk_ParquetAD4111!,
    (4, 2, 0, 1) => eval_sigmadk_ParquetAD4201!,
    (5, 0, 0, 1) => eval_sigmadk_ParquetAD5001!,
    (5, 0, 1, 1) => eval_sigmadk_ParquetAD5011!,
    (5, 1, 0, 1) => eval_sigmadk_ParquetAD5101!,
    (6, 0, 0, 1) => eval_sigmadk_ParquetAD6001!
)

end
