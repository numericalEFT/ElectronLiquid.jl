using ElectronLiquid
using FeynmanDiagram

diagtype = :sigma
order = 6
para = UEG.ParaMC(rs=1.0, beta=25, order=order, isDynamic=false)
_partition = UEG.partition(order)

partition = Vector{NTuple{4,Int}}()
for p in _partition
    push!(partition, (p..., 1))
end
FeynGraphs = Diagram.diagram_parquet_noresponse(:sigma, para, partition, [pr -> pr isa FrontEnds.BareGreenId,
        pr -> pr isa FrontEnds.BareInteractionId, pr -> pr.extK[1] != 0], optimize_level=1)

Diagram.compileC_ParquetAD_toFiles(FeynGraphs, order + 1, "sigmadk", compiler="icc")