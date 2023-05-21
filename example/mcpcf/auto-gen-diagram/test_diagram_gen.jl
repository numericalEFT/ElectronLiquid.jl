using ElectronLiquid
using FeynmanDiagram

const rs = 2.0
const beta = 25

const para = UEG.ParaMC(rs=rs, beta=beta)
const kF = para.kF

# partition = [(1, 0, 0), (2, 0, 0),]
partition = [(2, 0, 0),]
channel = [PHr,] # multiply by 2 for PHEr
# neighbor = UEG.neighbor(partition)

filter = [NoHartree, NoBubble]
# set proper to reduce ver3 diagram(Y).
# set transfer momentum for proper

diagram = Ver4.diagram(para, partition; channel=channel, filter=filter)
partition, diagpara, diag, root, extT = diagram

println(diagram)
println(diag)
println(diagpara)

const idx = 1

d = diag[idx]
dp = diagpara[idx]
Kn = dp.innerLoopNum
Tn = dp.totalTauNum - 1
println((Kn, Tn))
K = zeros(Float64, (3, 4))
K[1, 1], K[1, 2] = kF, kF
K[:, 3] = [kF, kF, kF] ./ sqrt(3)
T = zeros(Float64, 4)
weight = d.node.current
ExprTree.evalKT!(d, K, T, para)
w = sum(weight[r] for (ri, r) in enumerate(d.root))
println(w)
