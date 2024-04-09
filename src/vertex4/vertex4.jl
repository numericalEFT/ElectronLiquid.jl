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
include("ver4_lavg_beta.jl")
include("ver4_OAA.jl")
include("ver4_OAA_Clib.jl")

# include("ver4_PH_l_vegas.jl")
# include("ver4_PH_l_mcmc.jl")
# include("ver4_ParquetAD_compile_dynamic.jl")
# include("ver4KW.jl")

include("source_codeParquetAD/Cwrapper_vertex4_ParquetAD.jl")
include("source_codeParquetAD_Proper/Cwrapper_vertex4_ParquetAD.jl")

const evalfunc_vertex4_map = Dict(
    (0, 0, 0) => eval_vertex4_ParquetAD000!,
    (0, 0, 1) => eval_vertex4_ParquetAD001!,
    (0, 0, 2) => eval_vertex4_ParquetAD002!,
    (0, 0, 3) => eval_vertex4_ParquetAD003!,
    (0, 0, 4) => eval_vertex4_ParquetAD004!,
    (1, 0, 0) => eval_vertex4_ParquetAD100!,
    (1, 0, 1) => eval_vertex4_ParquetAD101!,
    (1, 0, 2) => eval_vertex4_ParquetAD102!,
    (1, 0, 3) => eval_vertex4_ParquetAD103!,
    (1, 1, 0) => eval_vertex4_ParquetAD110!,
    (1, 1, 1) => eval_vertex4_ParquetAD111!,
    (1, 1, 2) => eval_vertex4_ParquetAD112!,
    (1, 2, 0) => eval_vertex4_ParquetAD120!,
    (1, 2, 1) => eval_vertex4_ParquetAD121!,
    (1, 3, 0) => eval_vertex4_ParquetAD130!,
    (2, 0, 0) => eval_vertex4_ParquetAD200!,
    (2, 0, 1) => eval_vertex4_ParquetAD201!,
    (2, 0, 2) => eval_vertex4_ParquetAD202!,
    (2, 1, 0) => eval_vertex4_ParquetAD210!,
    (2, 1, 1) => eval_vertex4_ParquetAD211!,
    (2, 2, 0) => eval_vertex4_ParquetAD220!,
    (3, 0, 0) => eval_vertex4_ParquetAD300!,
    (3, 0, 1) => eval_vertex4_ParquetAD301!,
    (3, 1, 0) => eval_vertex4_ParquetAD310!,
    (4, 0, 0) => eval_vertex4_ParquetAD400!
)

const evalfunc_vertex4Proper_map = Dict(
    (0, 0, 0) => eval_vertex4Proper_ParquetAD000!,
    (0, 0, 1) => eval_vertex4Proper_ParquetAD001!,
    (0, 0, 2) => eval_vertex4Proper_ParquetAD002!,
    (0, 0, 3) => eval_vertex4Proper_ParquetAD003!,
    (0, 0, 4) => eval_vertex4Proper_ParquetAD004!,
    (0, 0, 5) => eval_vertex4Proper_ParquetAD005!,
    (1, 0, 0) => eval_vertex4Proper_ParquetAD100!,
    (1, 0, 1) => eval_vertex4Proper_ParquetAD101!,
    (1, 0, 2) => eval_vertex4Proper_ParquetAD102!,
    (1, 0, 3) => eval_vertex4Proper_ParquetAD103!,
    (1, 0, 4) => eval_vertex4Proper_ParquetAD104!,
    (1, 1, 0) => eval_vertex4Proper_ParquetAD110!,
    (1, 1, 1) => eval_vertex4Proper_ParquetAD111!,
    (1, 1, 2) => eval_vertex4Proper_ParquetAD112!,
    (1, 1, 3) => eval_vertex4Proper_ParquetAD113!,
    (1, 2, 0) => eval_vertex4Proper_ParquetAD120!,
    (1, 2, 1) => eval_vertex4Proper_ParquetAD121!,
    (1, 2, 2) => eval_vertex4Proper_ParquetAD121!,
    (1, 3, 0) => eval_vertex4Proper_ParquetAD130!,
    (1, 3, 1) => eval_vertex4Proper_ParquetAD131!,
    (1, 4, 0) => eval_vertex4Proper_ParquetAD140!,
    (2, 0, 0) => eval_vertex4Proper_ParquetAD200!,
    (2, 0, 1) => eval_vertex4Proper_ParquetAD201!,
    (2, 0, 2) => eval_vertex4Proper_ParquetAD202!,
    (2, 0, 3) => eval_vertex4Proper_ParquetAD203!,
    (2, 1, 0) => eval_vertex4Proper_ParquetAD210!,
    (2, 1, 1) => eval_vertex4Proper_ParquetAD211!,
    (2, 1, 2) => eval_vertex4Proper_ParquetAD212!,
    (2, 2, 0) => eval_vertex4Proper_ParquetAD220!,
    (2, 2, 1) => eval_vertex4Proper_ParquetAD221!,
    (2, 3, 0) => eval_vertex4Proper_ParquetAD230!,
    (3, 0, 0) => eval_vertex4Proper_ParquetAD300!,
    (3, 0, 1) => eval_vertex4Proper_ParquetAD301!,
    (3, 0, 2) => eval_vertex4Proper_ParquetAD302!,
    (3, 1, 0) => eval_vertex4Proper_ParquetAD310!,
    (3, 1, 1) => eval_vertex4Proper_ParquetAD311!,
    (3, 2, 0) => eval_vertex4Proper_ParquetAD320!,
    (4, 0, 0) => eval_vertex4Proper_ParquetAD400!,
    (4, 0, 1) => eval_vertex4Proper_ParquetAD401!,
    (4, 1, 0) => eval_vertex4Proper_ParquetAD410!,
    (5, 0, 0) => eval_vertex4Proper_ParquetAD500!
)

end
