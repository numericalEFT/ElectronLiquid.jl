using ElectronLiquid
using CompositeGrids
using FiniteDifferences
using JLD2

"""
    function KO(para::ParaMC, kamp=para.kF, kamp2=para.kF; N=100, mix=1.0, verbose=1, ct=false)
    
    Calculate Fs with the self-consistent equation
    ```math
    f_s = -\\left< (v_q+f_s)/(1-(v_q+f_s)Π_0) \\right>
    ```
"""
function KO(para::ParaMC, kamp=para.kF, kamp2=para.kF; N=100, mix=0.8, verbose=1, eps=1e-5)
    # Fp, Fm = para.Fs, para.Fa
    Fp, Fm = 0.0, 0.0
    for i = 1:N
        p_l = ParaMC(rs=para.rs, beta=para.beta, Fs=Fp, Fa=Fm, order=para.order, mass2=para.mass2, isDynamic=true, isFock=false)
        # println(p_l.Fs)
        wp, wm, angle = Ver4.exchange_interaction(p_l, kamp, kamp2; ct=false, verbose=verbose)
        Fp_new = -Ver4.Legrendre(0, wp, angle)
        if abs(Fp_new - Fp) < eps
            break
        else
            if verbose > 1
                println("iter $i: $Fp_new vs $Fp")
            end
        end
        Fp = mix * Fp_new + (1 - mix) * Fp
        Fm = -0.0
    end
    if verbose > 0
        println("Self-consistent approach: ")
        println("Fs = ", Fp)
        println("Fa = ", Fm)
    end
    return Fp, Fm
end

function dKO_dΛ(paras, Λgrid)
    # Fs = [KO(paras[li], lambda, lambda; verbose=0)[1] for lambda in Λgrid]
    dFs = -∂Rs_∂Λ_exchange(paras, Λgrid; ct=false) ./ (1 .+ ∂Rs_∂fs_exchange(paras, Λgrid; ct=false) ./ paras[1].NF)
    dFm = zero(Λgrid.grid)

    # F0 = -Interp.integrate1D(dFs, Λgrid, [Λgrid[1], Λgrid[end]])
    # Fs = KO(paras[1], Λgrid[1], Λgrid[1]; verbose=0, ct=ct)[1]
    # println(F0, " vs ", Fs[1])
    return dFs, dFm
end

"""
    function ∂Rs_∂Λ_exchange(paras, Λgrid; ct=false)
    
     return N_F ∂<R_{k_1-k_2}>/∂Λ
"""
function ∂Rs_∂Λ_exchange(paras, Λgrid; ct=false)

    function exchange(p, kamp, kamp2)
        wp, wm, angle = Ver4.exchange_interaction(p, kamp, kamp2; ct=ct)
        return Ver4.Legrendre(0, wp, angle)
    end

    dF = zero(Λgrid)
    for li in eachindex(Λgrid)
        lambda = Λgrid[li]
        p_l = paras[li]
        dFp = central_fdm(5, 1)(λ -> exchange(p_l, λ, λ), lambda) #use central finite difference method to calculate the derivative
        dF[li] = dFp
    end
    return dF
end

"""
    function ∂Rs_∂fs_exchange(paras, Λgrid; ct=false)
    
     return N_F ∂<R_{k_1-k_2}>/∂f
     
"""
function ∂Rs_∂fs_exchange(paras, Λgrid; ct=false)
    dF = zero(Λgrid)
    for li in eachindex(Λgrid)
        lambda = Λgrid[li]
        p_l = paras[li]
        wp, wm, angle = Ver4.exchange_interaction_df(p_l, lambda, lambda; ct=ct)
        dF[li] = Ver4.Legrendre(0, wp, angle)
    end
    return dF
end


function gamma4_treelevel_RG(para, Λgrid; verbose=1, rtol=1e-4, mix=0.9)

    fs = [KO(para, lambda, lambda; verbose=0)[1] for lambda in Λgrid]
    us = deepcopy(fs)

    idx = 1
    while true
        paras = [UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs[li], Fa=0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false) for li in eachindex(Λgrid)]

        ∂R_∂Λ = ∂Rs_∂Λ_exchange(paras, Λgrid; ct=false)
        ∂R_∂f = ∂Rs_∂fs_exchange(paras, Λgrid; ct=false) / para.NF

        dfs = -∂R_∂Λ ./ ∂R_∂f ./ 2.0
        dus = ∂R_∂Λ ./ 2.0
        fs_new = [-Interp.integrate1D(dfs, Λgrid, [Λgrid[idx], Λgrid[end]]) for idx in eachindex(Λgrid)]
        us_new = [-Interp.integrate1D(dus, Λgrid, [Λgrid[idx], Λgrid[end]]) for idx in eachindex(Λgrid)]

        max_fs = maximum(abs.((fs .- fs_new)))
        max_us = maximum(abs.((us .- us_new)))
        if verbose > 0
            println("iteration $(idx) with max_fs = $(max_fs) and max_us = $(max_us)")
        end
        if (max_fs / maximum(abs.(fs)) < rtol) && (max_us / maximum(abs.(us)) < rtol)
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

        @. fs = fs * (1 - mix) + fs_new * (mix)
        @. us = us * (1 - mix) + us_new * (mix)

        idx += 1
    end
    error("failed to converge!")
end

function gamma4_treelevel_RG_2(para, Λgrid; verbose=1, rtol=1e-4, mix=0.9)

    fs = [KO(para, lambda, lambda; verbose=0)[1] for lambda in Λgrid]
    us = deepcopy(fs)

    idx = 1
    while true
        paras = [UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs[li], Fa=0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false) for li in eachindex(Λgrid)]

        ∂R_∂Λ = ∂Rs_∂Λ_exchange(paras, Λgrid; ct=false)
        ∂R_∂f = ∂Rs_∂fs_exchange(paras, Λgrid; ct=false) / para.NF

        dfs = -∂R_∂Λ ./ 2.0
        dus = ∂R_∂Λ ./ 2.0
        fs_new = [-Interp.integrate1D(dfs, Λgrid, [Λgrid[idx], Λgrid[end]]) for idx in eachindex(Λgrid)]
        us_new = [-Interp.integrate1D(dus, Λgrid, [Λgrid[idx], Λgrid[end]]) for idx in eachindex(Λgrid)]

        max_fs = maximum(abs.((fs .- fs_new)))
        max_us = maximum(abs.((us .- us_new)))
        if verbose > 0
            println("iteration $(idx) with max_fs = $(max_fs) and max_us = $(max_us)")
        end
        if (max_fs / maximum(abs.(fs)) < rtol) && (max_us / maximum(abs.(us)) < rtol)
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

        @. fs = fs * (1 - mix) + fs_new * (mix)
        @. us = us * (1 - mix) + us_new * (mix)

        idx += 1
    end
    error("failed to converge!")
end


function gamma4_treelevel_KO(para, Λgrid; verbose=1)
    fs = [KO(para, lambda, lambda; verbose=0)[1] for lambda in Λgrid]
    paras = [UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs[li], Fa=0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false) for li in eachindex(Λgrid)]
    dfs = dKO_dΛ(paras, Λgrid)[1]
    us = -deepcopy(fs)
    dus = -deepcopy(us)
    if verbose > 0
        kF_idx = searchsortedfirst(Λgrid, para.kF)
        println("kF_idx: ", kF_idx, " with ", Λgrid[kF_idx] / para.kF, " kF")
        println("KO: Fs(kF) = -Us(kF) = ", fs[kF_idx])

        ######### test if dF_dΛ is correct #########
        F0 = -Interp.integrate1D(dfs, Λgrid, [Λgrid[1], Λgrid[end]])
        println("grid quality: $F0 vs $(fs[1])")
        @assert abs(F0 - fs[1]) < 1e-3 "dF_dΛ is not correct! $F0 vs $(fs[1])"
    end
    return fs, us, dfs, dus
end

if abspath(PROGRAM_FILE) == @__FILE__

    rs = [10.0,]
    mass2 = [1e-5,]
    beta = [25.0,]
    order = [1,]
    neval = 1e7

    for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)

        _para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=_order, mass2=_mass2, isDynamic=true)
        Λgrid = CompositeGrid.LogDensedGrid(:gauss, [1.0 * _para.kF, 100 * _para.kF], [_para.kF,], 8, 0.01 * _para.kF, 8)
        fs, us, dfs, dus = gamma4_treelevel_RG_2(_para, Λgrid; verbose=1)
        println(fs)
        println(dfs)

        println(us)
        println(dus)

        jldopen("data_f.jld2", "a+") do f
            key = "$(UEG.short(_para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (Λgrid, fs, us, dfs, dus)
        end
    end
end