using ElectronLiquid
using JLD2

dim = 3
rs = [1.0,]
mass2 = [0.5, 1.0, 2.0, 3.0]
Fs = [-0.0,]
beta = [25.0]
order = [3,]
neval = 1e6
isDynamic = false
isFock = false
spinPolarPara = 0.0 # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock)
    filename = "data_freeE.jld2"

    freeE, result = FreeEnergy.MC(para; neval=neval, filename=filename, spinPolarPara=spinPolarPara)
end