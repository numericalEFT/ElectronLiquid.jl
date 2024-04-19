using ElectronLiquid
using FeynmanDiagram

@inline function totalMomNum(order::Int, diagtype::Symbol)
    if diagtype == :vertex3
        return order + 2
    elseif diagtype == :vertex4
        return order + 3
    else
        return order + 1
    end
end

diagtype = :sigma # :sigma, :vertex3, :vertex4, :freeEnergy, :green, :chargePolar
order = 3
# filter = [Parquet.NoHartree, Parquet.Proper]
filter = [Parquet.NoHartree,]
KinL, KoutL, KinR = zeros(16), zeros(16), zeros(16)
KinL[1], KoutL[2], KinR[3] = 1.0, 1.0, 1.0

para = UEG.ParaMC(rs=1.0, beta=25, order=order, isDynamic=false)

if diagtype == :chargePolar || diagtype == :sigma
    _partition = UEG.partition(order)
else
    _partition = UEG.partition(order, offset=0)
end

if diagtype == :vertex4 || diagtype == :vertex3
    partition = Vector{NTuple{3,Int}}()
    for (o, sOrder, vOrder) in _partition
        o == 0 && sOrder > 0 && continue
        push!(partition, (o, sOrder, vOrder))
    end

    FeynGraphs = Diagram.diagram_parquet_response(diagtype, para, partition, optimize_level=1, filter=filter, transferLoop=KinL - KoutL)
elseif diagtype == :green || diagtype == :freeEnergy || diagtype == :chargePolar
    partition = Vector{NTuple{3,Int}}()
    for (o, sOrder, vOrder) in _partition
        o == 0 && vOrder > 0 && continue
        push!(partition, (o, sOrder, vOrder))
    end
    if diagtype == :freeEnergy
        FeynGraphs = Diagram.diagram_GV_freeE(para, partition, optimize_level=1)
    else
        FeynGraphs = Diagram.diagram_parquet_noresponse(diagtype, para, partition, optimize_level=1)
    end
else
    partition = _partition
    FeynGraphs = Diagram.diagram_parquet_noresponse(diagtype, para, partition, optimize_level=1)
end

# compile Python
for p in FeynGraphs[1]
    Compilers.compile_Python(FeynGraphs[3][p], "func_$(diagtype)_o$(p[1])$(p[2])$(p[3]).py")
end

println("Compile finished.")