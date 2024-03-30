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
import ..FeynmanDiagram.Parquet: DiagPara
using ..Measurements

# push!(LOAD_PATH, "../common/")
using ..UEG
using ..Propagator
# using ..Sigma
# import ..ExprTreeF64

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

function diagPara(para::ParaMC, order, filter, transferLoop)
    inter = [Parquet.Interaction(ChargeCharge, para.isDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    return DiagPara(
        type=Parquet.Ver4Diag,
        # innerLoopNum=order - 1,
        innerLoopNum=order,
        hasTau=true,
        spin=para.spin,
        firstLoopIdx=4,
        interaction=inter,
        filter=filter,
        transferLoop=transferLoop
    )
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

# include("common.jl")

# include("ver4_avg.jl")
# include("ver4KW.jl")
# include("ver4_PH_l.jl")
# include("ver4_PH_l_df.jl")
include("ver4_PH_l_AD.jl")
include("ver4_PH_l_AD_Clib.jl")
# include("ver4_generic.jl")
# include("exchange_interaction.jl")
# include("ver4_PP_l.jl")
# include("ver4_PH_l_vegas.jl")
include("ver4_PH_l_mcmc.jl")


include("ver4_ParquetAD.jl")
# include("ver4_ParquetAD_compile_dynamic.jl")
include("parquetAD_Clib.jl")

end