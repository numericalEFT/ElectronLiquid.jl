# include("../common/interaction.jl")

# using ElectronGas: Polarization
# using ElectronGas: SelfEnergy
# using GreenFunc
using ElectronGas
using Printf

beta = 25.0
mass2 = 0.001

# rslist = [5.0, 5.0, 5.0, 5.0, 5.0,]
rslist = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0]
# fplist = [-0.20633, -0.33343, -0.43340, -0.51587, -0.58545, -0.64494, -0.74100] #order 1, self-consistent
# fplist = [-0.22484, -0.38465, -0.52675, -0.65879, -0.78412, -0.90474, -1.1344] #order 1, variational outcome
# fplist = [-0.44, -0.76, -1.0, -1.2, -1.48, -1.72, -2.2] #varitional parameter
fplist = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

paralist = [Parameter.rydbergUnit(1.0 / beta, rs, 3, Λs=mass2) for rs in rslist]

for (idx, rs) in enumerate(rslist)
    fp = fplist[idx]
    para = paralist[idx]

    Σ = SelfEnergy.G0W0(para; Euv=100 * para.EF, maxK=8 * para.kF, Nk=16, order=8, minK=1e-8 * para.kF, int_type=:ko_const, Fs=-fp, Fa=-0.0)
    zz = SelfEnergy.zfactor(Σ)
    dS_dw = 1 - 1 / zz

    # println("zfactor ")
    # println(zz)
    # println(dS_dw)

    # effective mass is not stable!
    # if you use two different defintion of k1 and k2, one got different effective mass
    # println("mass ratio")
    ratio = SelfEnergy.massratio(para, Σ)
    # println(ratio)
    dS_dK = (1 / ratio / zz - 1)
    # println(dS_dK)
    dmu = SelfEnergy.chemicalpotential(para, Σ)
    @printf("%4i    %16.8f  %16.8f  %16.8f  %16.8f  %16.8f\n", rs, zz, ratio, dS_dw, dS_dK, dmu)
end

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