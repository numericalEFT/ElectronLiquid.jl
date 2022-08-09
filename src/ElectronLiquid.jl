module ElectronLiquid

using StaticArrays
using Measurements
using MCIntegration
using FeynmanDiagram
export TwoBodyChannel, Alli, PHr, PHEr, PPr, AnyChan

using CompositeGrids
using Lehmann
using ElectronGas
using Parameters

const ExprTreeF64 = ExpressionTree{ExprTree.LoopPool{Float64},DiagramId,Float64,Float64}

include("./common/para_builder.jl")
using .UEG
export UEG

include("./common/counterterm.jl")
using .CounterTerm
export CounterTerm

include("./common/eval.jl")
using .Propagator

include("./sigma/sigma.jl")
using .Sigma
export Sigma

include("./vertex4/vertex4.jl")
using .Ver4
export Ver4


end