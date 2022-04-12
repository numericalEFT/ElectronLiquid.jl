include("../common/interaction.jl")

using ElectronGas: Polarization
using ElectronGas: SelfEnergy
using GreenFunc

Σ = SelfEnergy.G0W0(para; Euv=100 * para.EF, maxK=8 * para.kF, Nk=12, order=6, minK=1e-8 * para.kF, int_type=:ko_const, Fs=Fs)
zz = SelfEnergy.zfactor(Σ)
println(zz)

# effective mass is not stable!
# if you use two different defintion of k1 and k2, one got different effective mass
println(SelfEnergy.massratio(para, Σ))

# kgrid = Σ.spaceGrid
# kF_label = searchsortedfirst(kgrid.grid, kF)
# Σ_freq = GreenFunc.toMatFreq(Σ, [0, 1])
# k1 = kF_label
# println("kF = ", para.kF)
# println("k0 = ", kgrid.grid[k1])

# for dk in 1:30
#     k2 = kF_label + dk
#     sigma1 = real(Σ_freq.dynamic[1, 1, k1, 1] + Σ_freq.instant[1, 1, k1])
#     sigma2 = real(Σ_freq.dynamic[1, 1, k2, 1] + Σ_freq.instant[1, 1, k2])
#     ds_dk = (sigma1 - sigma2) / (kgrid.grid[k1] - kgrid.grid[k2])

#     # println("m/kF ds_dk = $(me/kF*ds_dk)")
#     println(kgrid.grid[k2] - kgrid.grid[k1], "    ", 1.0 / zz / (1 + me / kF * ds_dk))
# end