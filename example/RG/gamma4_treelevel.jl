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

function gamma4_treelevel_RG(para, Λgrid; verbose=1, rtol=1e-4, mix=0.9)
    fs = [KO(para, lambda, lambda; verbose=0)[1] for lambda in Λgrid]
    us = deepcopy(fs)

    dfs = zero(Λgrid.grid)
    dus = zero(Λgrid.grid)

    idx = 1
    while true
        fs_new, us_new = zero(fs), zero(us)
        flag = true
        Threads.@threads for li in eachindex(Λgrid)
            lambda = Λgrid[li]
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
            # if abs((fs[li] - fs_new[li]) / fs[li]) > rtol || abs((us[li] - us_new[li]) / us[li]) > rtol
            #     flag = false
            # end
        end
        max_fs = maximum(abs.((fs .- fs_new)))
        max_us = maximum(abs.((us .- us_new)))
        if max_fs / maximum(fs) < rtol || max_us / maximum(us) < rtol
            flag = break
        end

        if verbose > 0
            println("iteration $(idx) with max_fs = $(max_fs) and max_us = $(max_us)")
        end

        @. fs = fs * (1 - mix) + fs_new * (mix)
        @. us = us * (1 - mix) + us_new * (mix)
        idx += 1
    end
    if verbose >= 0
        println("total iteration $(idx)")
    end
    if verbose > 0
        kF_idx = searchsortedfirst(Λgrid, para.kF)
        println("kF_idx: ", kF_idx, " with ", Λgrid[kF_idx] / para.kF, " kF")
        println("Fs(kF): ", fs[kF_idx])
        println("Us(kF): ", us[kF_idx])
    end
    return fs, us
end


function gamma4_treelevel_KO(para, Λgrid; verbose=1)
    fs = [-KO(para, lambda, lambda; verbose=0)[1] for lambda in Λgrid]
    us = -deepcopy(fs)
    if verbose > 0
        kF_idx = searchsortedfirst(Λgrid, para.kF)
        println("kF_idx: ", kF_idx, " with ", Λgrid[kF_idx] / para.kF, " kF")
        println("KO: Fs(kF) = -Us(kF) = ", fs[kF_idx])
    end
    return fs, us
end
