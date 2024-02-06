using ElectronLiquid, FeynmanDiagram
import ..FeynmanDiagram.FrontEnds: TwoBodyChannel, Alli, PHr, PHEr, PPr, AnyChan
import ..FeynmanDiagram.FrontEnds: Filter, NoHartree, NoFock, DirectOnly, Wirreducible, Girreducible, NoBubble, Proper
import ..FeynmanDiagram.FrontEnds: Response, Composite, ChargeCharge, SpinSpin, UpUp, UpDown
import ..FeynmanDiagram.FrontEnds: AnalyticProperty, Instant, Dynamic


const order = 4
const isDynamic = false

if !isDynamic
    para = UEG.ParaMC(rs=1.0, beta=25, order=order, isDynamic=false)
    _partition = UEG.partition(order, offset=0)
    partition = Vector{Tuple{Int64,Int64,Int64}}()
    for (o, sOrder, vOrder) in _partition
        o == 0 && sOrder > 0 && continue
        push!(partition, (o, sOrder, vOrder))
    end

    println("generating diagrams")
    diagram = Ver4.diagramParquet(para, partition; filter=[NoHartree,])
    println("diagram generated")
    partition, diagpara, FeynGraphs, extT_labels, spin_conventions = diagram

    MaxLoopNum = maximum([key[1] for key in partition]) + 3
    Ver4.compileC_ParquetAD_toFiles(para.order, partition, FeynGraphs, MaxLoopNum, compiler="icc")
else
    para = UEG.ParaMC(rs=1.0, beta=25, order=order, isDynamic=true)
    partition = UEG.partition(para.order)

    println("generating diagrams")
    diagram = Ver4.diagramParquet(para, partition; filter=[NoHartree, NoBubble])
    println("diagram generated")
    partition, diagpara, FeynGraphs, extT_labels, spin_conventions = diagram

    MaxLoopNum = maximum([key[1] for key in partition]) + 3
    Ver4.compileC_ParquetAD_toFiles_dynamic(para.order, partition, FeynGraphs, MaxLoopNum, compiler="icc")
end