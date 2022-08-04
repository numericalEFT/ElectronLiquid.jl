using ElectronLiquid
using FeynmanDiagram
p = (2,0,0)
mass2=0.01
para = ElectronLiquid.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=1, mass2=mass2, isDynamic=true)
diagram = Ver4.diagram(para, [p,]; filter=[Proper, NoBubble], channel=[PHr, PHEr, PPr])
data, result = Ver4.PH(para, diagram; neval=1e8, print=0, l=[0,], n=[-1, 0, 0])
if isnothing(data) == false
    println(data)
end