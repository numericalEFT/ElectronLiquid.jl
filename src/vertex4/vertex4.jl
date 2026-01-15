module Ver4

using JLD2, CSV
using Printf, LinearAlgebra, DataFrames
using ..StaticArrays
using ..Parameters
using ..CompositeGrids
using ..ElectronGas
using ..MCIntegration
using ..Lehmann

using ..FeynmanDiagram
import ..FeynmanDiagram.FrontEnds: TwoBodyChannel, Alli, PHr, PHEr, PPr, AnyChan
import ..FeynmanDiagram.FrontEnds: Filter, NoHartree, NoFock, DirectOnly, Wirreducible, Girreducible, NoBubble, Proper
import ..FeynmanDiagram.FrontEnds: Response, Composite, ChargeCharge, SpinSpin, UpUp, UpDown
import ..FeynmanDiagram.FrontEnds: AnalyticProperty, Instant, Dynamic
import ..FeynmanDiagram.Parquet: DiagPara, Ver4Diag
using ..Measurements
using ..Diagram
# push!(LOAD_PATH, "../common/")
using ..UEG
using ..Propagator
import ..Propagator: LeafStateADDynamic

import ..Weight

"""
    struct OneAngleAveraged

The parameters for the one-angle-averaged vertex4.

# Members
- `para`: the parameters for the MC integration
- `kamp`: the amplitude of the external momentum: [left_leg, right_legs]
- `ωn`: vector of the frequency of the external legs, each element is a 3-vector [left_in, left_out, right_in]
- `channel`: the channel of the vertex4, :PH or :PP
- `l`: the angular momentum of the angle average
"""
struct OneAngleAveraged
    para::ParaMC
    kamp::Vector{Float64}
    ωn::Vector{Vector{Int}} #allow measure multiple frequency simultaneously
    channel::Symbol #:PH or :PP
    l::Int #angular momentum
    function OneAngleAveraged(para, kamp, ωn, channel, l)
        @assert channel == :PH || channel == :PP "the channel should be :PH or :PP"
        @assert length(kamp) == 2 "there two amplitude of K"
        # @assert length(ωn) == 3 "the length of ωn should be 3, which corresponds to Lin, Rin, Lout."
        return new(para, kamp, ωn, channel, l)
    end
end

function diagram_loadinfo(paramc::ParaMC, _partition::Vector{T};
    filter=[NoHartree], transferLoop=nothing,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"), filename="extvars_vertex4.jld2"
) where {T}
    diagpara = Vector{DiagPara}()
    extT_labels = Vector{Vector{Int}}[]
    spin_conventions = Vector{Response}[]

    fname = joinpath(root_dir, filename)
    jldopen(fname, "r") do f
        for p in _partition
            key_str = join(string.(p))
            if key_str in keys(f)
                extT, ext_spin = f[key_str]
                push!(diagpara, Diagram.diagPara(Ver4Diag, paramc.isDynamic, p[1], paramc.spin, filter, transferLoop))
                push!(extT_labels, extT)
                push!(spin_conventions, ext_spin)
            else
                error("$(key_str) not found in $(fname)")
            end
        end
    end
    return (_partition, diagpara, extT_labels, spin_conventions)
end

@inline function legendfactor(x, l, dim)
    if dim == 3
        if l == 0
            factor = 0.5
        elseif l == 1
            factor = x / 2.0
        elseif l == 2
            factor = (3x^2 - 1) / 4.0
        elseif l == 3
            factor = (5x^3 - 3x) / 4.0
        elseif l == 4
            factor = (35x^4 - 30x^2 + 3) / 16.0
        elseif l == 5
            factor = (63x^5 - 70x^3 + 15x) / 16.0
        else
            error("not implemented for $l channel in $dim-D")
        end
    elseif dim == 2
        factor = cos(l * x) / 2π
    else
        error("not implemented in $dim-D")
    end
    return factor
end

@inline function phase(varT, extT, ninL, noutL, ninR, β)
    # println(extT)
    tInL, tOutL, tInR, tOutR = varT[extT[INL]], varT[extT[OUTL]], varT[extT[INR]], varT[extT[OUTR]]
    winL, woutL, winR = π * (2ninL + 1) / β, π * (2noutL + 1) / β, π * (2ninR + 1) / β
    woutR = winL + winR - woutL
    return exp(-1im * (tInL * winL - tOutL * woutL + tInR * winR - tOutR * woutR))
end

@inline function phase(varT, extT, n, β)
    # println(extT)
    return phase(varT, extT, n[1], n[2], n[3], β)
end

@inline ud2sa(Wuu, Wud) = @. (Wuu + Wud) / 2, (Wuu - Wud) / 2
@inline sa2ud(Ws, Wa) = @. Ws + Wa, Ws - Wa

include("exchange_interaction.jl")
include("ver4_lavg.jl")
include("ver4_lavg_Clib.jl")
include("ver4_lavg_Project.jl")
include("ver4_lavg_beta.jl")
include("ver4_OAA.jl")
include("ver4_OAA_Clib.jl")
include("ver4_AR.jl")
include("ver4_Spec.jl")
include("ver4_Spec_Jl.jl")
# include("ver4_PH_l_vegas.jl")
# include("ver4_PH_l_mcmc.jl")
# include("ver4_ParquetAD_compile_dynamic.jl")
# include("ver4KW.jl")

include("source_codeParquetAD/Cwrapper_vertex4_ParquetAD.jl")
# include("source_codeParquetAD_Proper/Cwrapper_vertex4_ParquetAD.jl")
include("source_codeParquetAD_Proper_NoAlli/Cwrapper_vertex4_ParquetAD.jl")

const evalfunc_vertex4_map = Dict{Tuple{Int,Int,Int},Function}()
const evalfunc_vertex4Proper_map = Dict{Tuple{Int,Int,Int},Function}()

for sym in names(@__MODULE__; all=true)
    s_str = String(sym)

    m_ver4 = match(r"^eval_vertex4_ParquetAD(\d)(\d)(\d)!$", s_str)
    if m_ver4 !== nothing
        o, i, j = parse.(Int, m_ver4.captures)
        evalfunc_vertex4_map[(o, i, j)] = getfield(@__MODULE__, sym)
        continue
    end

    m_ver4Proper = match(r"^eval_vertex4Proper_ParquetAD(\d)(\d)(\d)!$", s_str)
    if m_ver4Proper !== nothing
        o, i, j = parse.(Int, m_ver4Proper.captures)
        evalfunc_vertex4Proper_map[(o, i, j)] = getfield(@__MODULE__, sym)
    end
end

end
