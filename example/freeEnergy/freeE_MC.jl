using ElectronLiquid
using JLD2

dim = 2
rs = [1.0,]
mass2 = [2.0,]
Fs = [-0.0,]
beta = [25.0]
order = [3,]
neval = 1e6
isDynamic = false
isFock = false

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock)
    filename = "data_freeE.jld2"

    freeE, result = FreeEnergy.MC(para; neval=neval, filename=filename)
end