using ElectronLiquid
using ElectronGas
using CompositeGrids
using FiniteDifferences
using GreenFunc

function zfactor_oneloop_RG(para, Λgrid, fs; verbose=1, rtol=1e-3, mix=1.0, ngrid=[0, 1])
    dz = zero(Λgrid.grid)
    z = zero(Λgrid.grid)
    dk = para.kF * 1e-2

    Threads.@threads for li in eachindex(Λgrid[1:end-1])
        lambda = Λgrid[li]
        p_l = UEG.ParaMC(rs=rs, beta=beta, Fs=fs[li], Fa=0.0, order=1, mass2=mass2, isDynamic=true, isFock=false)
        sigmadyn1, sigmainst1 = SelfEnergy.G0W0(p_l.basic, [lambda,]; maxK=Λgrid[end])
        sigmadyn2, sigmainst2 = SelfEnergy.G0W0(p_l.basic, [lambda + dk,]; maxK=Λgrid[end])

        sigmadyn1 = dlr_to_imfreq(to_dlr(sigmadyn1), ngrid)
        sigmadyn2 = dlr_to_imfreq(to_dlr(sigmadyn2), ngrid)
        zk1 = imag(sigmadyn1[1, 1] - sigmadyn1[2, 1]) / (2 * π * (ngrid[1] - ngrid[2]) / para.basic.β)
        zk2 = imag(sigmadyn2[1, 1] - sigmadyn2[2, 1]) / (2 * π * (ngrid[1] - ngrid[2]) / para.basic.β)
        println(lambda, " ", zk1, " ", zk2)
        dz[li] = (zk2 - zk1) / dk
    end

    println(dz)
    for li in eachindex(Λgrid)
        z[li] = Interp.integrate1D(dz, Λgrid, [Λgrid[li], Λgrid[end]])
    end

    if verbose > 0
        kF_idx = searchsortedfirst(Λgrid, para.kF)
        println("kF_idx: ", kF_idx, " with ", Λgrid[kF_idx] / para.kF, " kF")
        println("z(kF): ", z[kF_idx])
    end
    return z
end