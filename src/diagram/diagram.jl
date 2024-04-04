module Diagram

using JLD2, CSV
using Printf, LinearAlgebra, DataFrames
using ..StaticArrays
using ..Parameters

using ..FeynmanDiagram
import ..FeynmanDiagram.FrontEnds: TwoBodyChannel, Alli, PHr, PHEr, PPr, AnyChan
import ..FeynmanDiagram.FrontEnds: Filter, NoHartree, NoFock, DirectOnly, Wirreducible, Girreducible, NoBubble, Proper
import ..FeynmanDiagram.FrontEnds: Response, Composite, ChargeCharge, SpinSpin, UpUp, UpDown
import ..FeynmanDiagram.FrontEnds: AnalyticProperty, Instant, Dynamic
import ..FeynmanDiagram.Parquet: DiagPara, Interaction
import ..FeynmanDiagram.Parquet: DiagramType, VacuumDiag, SigmaDiag, GreenDiag, PolarDiag, Ver3Diag, Ver4Diag

using ..UEG

include("generate.jl")
include("compile.jl")

end