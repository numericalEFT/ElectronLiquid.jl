module ElectronLiquid

using StaticArrays
using Measurements
using MCIntegration
export report

using FeynmanDiagram
export TwoBodyChannel, Alli, PHr, PHEr, PPr, AnyChan

using CompositeGrids
using Lehmann
using ElectronGas
using Parameters

mutable struct Weight{T} <: FieldVector{2,T}
    d::T
    e::T
    Weight{T}() where {T} = new{T}(0.0, 0.0)
    Weight(d::T, e::T) where {T} = new{T}(d, e)
end

const Base.zero(::Type{Weight}) = Weight(0.0, 0.0)
const Base.abs(w::Weight) = abs(w.d) + abs(w.e) # define abs(Weight)
export Weight

const ExprTreeF64 = FeynmanDiagram.ExpressionTree{FeynmanDiagram.ExprTree.LoopPool{Float64},DiagramId,Float64,Float64}

include("./common/para_builder.jl")
using .UEG
export UEG
export ParaMC, getK

include("./common/renormalization.jl")
using .Renorm
export Renorm

include("./common/counterterm.jl")
using .CounterTerm
export CounterTerm

include("./common/eval.jl")
using .Propagator

include("./green/green.jl")
using .Green
export Green

include("./sigma/sigma.jl")
using .Sigma
export Sigma

include("./vertex4/vertex4.jl")
using .Ver4
export Ver4

include("./vertex3/vertex3.jl")
using .Ver3
export Ver3

using SnoopPrecompile
@precompile_all_calls begin
    # In here put "toy workloads" that exercise the code you want to precompile
    # para = UEG.ParaMC(rs=5.0, beta=25.0)
    # Sigma.diagram(para, [(2, 0, 0),])
    # Ver4.diagram(para, [(2, 0, 0),])
end

end