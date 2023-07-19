using Test
using ElectronLiquid
using FeynmanDiagram
using FiniteDifferences
using Lehmann
using Measurements

function compare(data, expect)
    # println(data, ", ", expect)
    @test isapprox(data.val, expect, atol=5 * data.err)
end

if isempty(ARGS)
    include("para.jl")
    include("interaction.jl")
    include("renormalization.jl")
    include("counterterm.jl")
    include("vertex4.jl")
    include("vertex3.jl")
    include("RG.jl")
else
    include(ARGS[1])
end
