using ElectronLiquid
using CompositeGrids
using FiniteDifferences
using JLD2

function Fs(para, kamp, kamp2; ct=false)
    # wp, wm, angle = Ver4.exchange_interaction(para, kamp, kamp2; ct=ct)
    Ws0, Wa0 = Ver4.projected_exchange_interaction(0, para, Ver4.exchange_interaction; kamp=kamp, kamp2=kamp2, verbose=0, ct=ct)
    return -Ws0
    # return -Ver4.Legrendre(0, wp, angle)
end

function Fa(para, kamp, kamp2; ct=false)
    Ws0, Wa0 = Ver4.projected_exchange_interaction(0, para, Ver4.exchange_interaction; kamp=kamp, kamp2=kamp2, verbose=0, ct=ct)
    return -Wa0
end

"""
    function KO(para::ParaMC, kamp=para.kF, kamp2=para.kF; N=100, mix=1.0, verbose=1, ct=false)
    
    Calculate Fs with the self-consistent equation
    ```math
    f_s = -\\left< (v_q+f_s)/(1-(v_q+f_s)Π_0) \\right>
    ```
"""
function KO(para::ParaMC, kamp=para.kF, kamp2=para.kF; N=100, mix=1.0, verbose=1, ct=false)
    Fp, Fm = para.Fs, para.Fa
    for i = 1:N
        p_l = ParaMC(rs=para.rs, beta=para.beta, Fs=Fp, Fa=Fm, order=para.order, mass2=para.mass2, isDynamic=true, isFock=false)
        # println(p_l.Fs)
        wp, wm, angle = Ver4.exchange_interaction(p_l, kamp, kamp2; ct=false)
        Fp = -Ver4.Legrendre(0, wp, angle)
        Fm = -0.0
    end
    if verbose > 0
        println("Self-consistent approach: ")
        println("Fs = ", Fp)
        println("Fa = ", Fm)
    end
    return Fp, Fm
end

function ∂Rs_∂Λ_exchange(paras, Λgrid; ct=false)
    dFs, dFa = zero(Λgrid.grid)
    for li in eachindex(Λgrid)
        lambda = Λgrid[li]
        p_l = paras[p_l]
        dFp = central_fdm(5, 1)(λ -> Fs(p_l, λ, λ; ct=ct), lambda) #use central finite difference method to calculate the derivative
        dFm = central_fdm(5, 1)(λ -> Fa(p_l, λ, λ; ct=ct), lambda) #use central finite difference method to calculate the derivative
        dFs[li] = dFp
        dFa[li] = dFm
    end
    return dFs, dFa
end

function ∂Rs_∂fs_exchange(paras, Λgrid; ct=false)
    dFs, dFa = zero(Λgrid.grid)
    for li in eachindex(Λgrid)
        lambda = Λgrid[li]
        p_l = paras[p_l]
        Ws0, Wa0 = Ver4.projected_exchange_interaction(0, p_l, Ver4.exchange_interaction_df; kamp=lambda, kamp2=lambda, verbose=0, ct=ct)
        dFs[li] = -Ws0
        dFa[li] = -Wa0
    end
    return dFs, dFa
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
        paras = [UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs[li], Fa=0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false) for li in eachindex(Λgrid)]
        for li in eachindex(Λgrid)
            lambda = Λgrid[li]
            p_l = paras[li]

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

    return fs, us, dfs, dus
end


function gamma4_treelevel_KO(para, Λgrid; verbose=1)
    fs = [KO(para, lambda, lambda; verbose=0)[1] for lambda in Λgrid]
    us = -deepcopy(fs)
    if verbose > 0
        kF_idx = searchsortedfirst(Λgrid, para.kF)
        println("kF_idx: ", kF_idx, " with ", Λgrid[kF_idx] / para.kF, " kF")
        println("KO: Fs(kF) = -Us(kF) = ", fs[kF_idx])
    end
    return fs, us
end

if abspath(PROGRAM_FILE) == @__FILE__

    rs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
    mass2 = [1e-5,]
    beta = [100.0,]
    order = [1,]
    neval = 1e7

    for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)

        _para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=_order, mass2=_mass2, isDynamic=true)
        Λgrid = CompositeGrid.LogDensedGrid(:gauss, [1.0 * para.kF, 20 * para.kF], [para.kF,], 8, 0.01 * para.kF, 8)
        fs, us, dfs, dus = gamma4_treelevel_RG(_para, Λgrid; verbose=1)

        jldopen("data_f.jld2", "a+") do f
            key = "$(UEG.short(_para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (_para, Λgrid, fs, us, dfs, dus)
        end
    end
end