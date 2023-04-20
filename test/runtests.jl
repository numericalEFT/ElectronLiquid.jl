using Test
using ElectronLiquid
using FeynmanDiagram
using FiniteDifferences
using Lehmann

if isempty(ARGS)
    include("interaction.jl")
    include("renormalization.jl")
    include("counterterm.jl")
    include("vertex4.jl")
    include("vertex3.jl")
else
    include(ARGS[1])
end
