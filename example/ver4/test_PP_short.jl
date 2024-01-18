using ElectronLiquid
order = 4
neval = 1e6
para = UEG.ParaMC(rs=1.0, beta=25, mass2=3.5, order=order, isDynamic=false)
# partition = [(o, 0, 0) for o in 1:order]
partition = UEG.partition(para.order)
ver4, result = Ver4.MC_PP_ParquetAD_Clib(para; partition=partition, neval=neval)
# ver4, result = Ver4.MC_PP_ParquetAD(para; partition=partition, neval=neval)
println(ver4)
