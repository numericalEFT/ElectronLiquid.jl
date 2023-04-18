using ElectronLiquid
using CompositeGrids
using FiniteDifferences

function Fs(para, kamp, kamp2)
    wp, wm, angle = Ver4.exchange_interaction(para, kamp, kamp2; ct=false)
    return -Ver4.Legrendre(0, wp, angle)
end

function Fa(para, kamp, kamp2)
    wp, wm, angle = Ver4.exchange_interaction(para, kamp, kamp2; ct=false)
    return -Ver4.Legrendre(0, wm, angle)
end

function KO(para::ParaMC, kamp=para.kF, kamp2=para.kF; N=100, mix=1.0, verbose=1)
    Fp, Fm = -para.Fs, -para.Fa
    for i = 1:N
        p_l = ParaMC(rs=para.rs, beta=para.beta, Fs=-Fp, Fa=-Fm, order=para.order, mass2=para.mass2, isDynamic=false, isFock=false)
        # println(p_l.Fs)
        wp, wm, angle = Ver4.exchange_interaction(p_l, kamp, kamp2; ct=false)
        nFs = Ver4.Legrendre(0, wp, angle)
        nFa = Ver4.Legrendre(0, wm, angle)
        # println("Fp = $nFs, Fold = $(p_l.Fs)")
        # Fp = Fp * (1 - mix) + nFs * mix # 
        # Fm = Fm * (1 - mix) + nFa * mix
        Fm = 0.0
        Fp = nFs
    end
    if verbose > 0
        println("Self-consistent approach: ")
        println("Fs = ", Fp)
        println("Fa = ", Fm)
    end
    return Fp, Fm
end

function gamma4_treelevel_RG(para, Λgrid)
    fs = [KO(para, lambda, lambda; verbose=0)[1] for lambda in Λgrid]
    us = deepcopy(fs)

    dfs = zero(Λgrid.grid)
    dus = zero(Λgrid.grid)

    idx = 1
    while true
        fs_new, us_new = zero(fs), zero(us)
        flag = true
        for (li, lambda) in enumerate(Λgrid)
            p_l = UEG.ParaMC(rs=rs, beta=beta, Fs=fs[li], Fa=0.0, order=1, mass2=mass2, isDynamic=true, isFock=false)

            # _Fs_dΛ = Fs(p_l, lambda + dΛ, lambda + dΛ)
            # _Fs = Fs(p_l, lambda, lambda)
            # dFs = (_Fs_dΛ - _Fs) / dΛ
            dFs = central_fdm(5, 1)(λ -> Fs(p_l, λ, λ), lambda) #use central finite difference method to calculate the derivative

            wp, wm, angle = Ver4.exchange_interaction_df(p_l, lambda, lambda; ct=false)
            Fs_df = Ver4.Legrendre(0, wp, angle) / para.NF

            dfs[li] = -dFs / Fs_df / 2.0
            dus[li] = dFs / 2.0

            fs_new[li] = Interp.integrate1D(dfs, Λgrid, [Λgrid[li], Λgrid[end]])
            us_new[li] = Interp.integrate1D(dus, Λgrid, [Λgrid[li], Λgrid[end]])
            if abs(fs[li] - fs_new[li]) > 1e-4 || abs(us[li] - us_new[li]) > 1e-4
                flag = false
            end
        end

        if flag
            break
        end
        fs .= fs_new
        us .= us_new
        idx += 1
    end
    println("iteration $(idx)")
    # kF_idx = 1
    # println("kF_idx: ", kF_idx, " with ", Λgrid[kF_idx] / para.kF, " kF")
    # println("Fs(kF): ", fs[kF_idx])
    # println("Us(kF): ", us[kF_idx])
    return fs, us
end


function gamma4_treelevel_KO(para, Λgrid)
    fs = [-KO(para, lambda, lambda; verbose=0)[1] for lambda in Λgrid]
    us = [KO(para, lambda, lambda; verbose=0)[1] for lambda in Λgrid]
    return fs, us
end
