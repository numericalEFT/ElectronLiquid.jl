using ElectronLiquid

para = ParaMC(rs=5.0, beta=25.0, mass2=0.001,
    Fs=-0.0, Fa=-0.0, isDynamic=true, order=2)
p = Ver4.OneAngleAveraged(para,
    [para.kF, para.kF], [[-1, 0, 0, -1],], :PP, 0)
diag = Ver4.diagram(para, UEG.partition(2);
    channel=[PPr, PHEr, PHr])
data, res = ElectronLiquid.Ver4.one_angle_averaged([p,], diag)
println(data[(2, 0, 0)])