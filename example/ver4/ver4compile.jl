using ElectronLiquid, FeynmanDiagram

para = UEG.ParaMC(rs=1.0, beta=25, order=6, isDynamic=false)
partition = UEG.partition(para.order)

println("generating diagrams")
diagram = Ver4.diagramParquet(para, partition; channel=channel = [PHr, PHEr, PPr,], filter=[NoHartree,])
println("diagram generated")
partition, diagpara, FeynGraphs, extT_labels, spin_conventions = diagram

MaxLoopNum = maximum([key[1] for key in partition]) + 2
Ver4.compileC_ParquetAD_toFiles(para.order, partition, FeynGraphs, MaxLoopNum)