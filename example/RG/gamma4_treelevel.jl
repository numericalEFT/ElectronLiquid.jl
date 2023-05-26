using ElectronLiquid
using CompositeGrids
using FiniteDifferences
using Roots
using JLD2

"""
    function KO(para::ParaMC, kamp=para.kF, kamp2=para.kF; N=100, mix=1.0, verbose=1, ct=false)
    
    Calculate Fs with the self-consistent equation
    ```math
    f_s = -\\left< (v_q+f_s)/(1-(v_q+f_s)Π_0) \\right>
    ```
"""
function KO(para::ParaMC, kamp=para.kF, kamp2=para.kF; a_s=0.0, N=100, mix=0.8, verbose=0, eps=1e-5, spin_spin=false, Fp=-0.2 * para.rs + 4π / para.me * a_s, Fm=0.0)

    u = 0.0

    function _u(fs)
        p_l = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs, Fa=0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false)
        wp, wm, angle = Ver4.exchange_interaction(p_l, kamp, kamp2; ct=false, verbose=verbose)
        return Ver4.Legrendre(0, wp, angle) + fs - 4 * π / para.me * a_s
    end

    function _u_f(fs)
        p_l = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs, Fa=0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false)
        wp, wm, angle = Ver4.exchange_interaction_df(p_l, kamp, kamp2; ct=false, verbose=verbose)
        return (Ver4.Legrendre(0, wp, angle)) / para.NF + 1.0
    end

    function newtown(fs)
        iter = 1
        err = eps * 10
        while err > eps && iter < N
            if verbose > 1
                println("$fs ->", fs - _u(fs) / _u_f(fs), " with ", _u(fs), ", ", _u_f(fs))
            end
            fx = _u(fs)
            fs_new = fs - fx / _u_f(fs)
            # err = abs(fs - fs_new)
            err = abs(fx)
            fs = fs_new
            iter += 1
            if iter >= N
                @warn("Newton-Raphson method doesn't converge. error = $err, got $fs")
            end
        end
        return fs
    end

    if spin_spin == false
        _Fs = newtown(Fp)
        p_l = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=_Fs, Fa=0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false)
        wp, wm, angle = Ver4.exchange_interaction(p_l, kamp, kamp2; ct=false, verbose=verbose)
        u = Ver4.Legrendre(0, wp, angle)
        if verbose > 0
            println("Self-consistent approach: ")
            println("Fs = ", _Fs)
            println("Fa = ", 0.0)
            println("u = ", u)
            println("4π a_s/m = ", 4π / para.me * a_s)
            # @assert abs(u + _Fs - 4π / para.me * a_s) < 10 * eps "u is not consistent with 4π a_s/m with error $(u + _Fs - 4π / para.me * a_s)"
        end
        return _Fs, 0.0, u
    else
        @assert spin_spin == false "Spin-Spin KO interaciton doesn't work yet."
    end

    # Fp, Fm = 0.0, 0.0
    # Fp, Fm = Fs, Fa
    # for i = 1:N
    #     p_l = ParaMC(rs=para.rs, beta=para.beta, Fs=Fp, Fa=Fm, order=para.order, mass2=para.mass2, isDynamic=true, isFock=false)
    #     # println(p_l.Fs)
    #     wp, wm, angle = Ver4.exchange_interaction(p_l, kamp, kamp2; ct=false, verbose=verbose)
    #     if spin_spin
    #         # u = 2<R^-> = -2<R^+>
    #         u = -Fp + 3Fm
    #         Fp_new = ufp(p_l, -u / 2, kamp, kamp2; verbose=verbose)
    #         Fm_new = ufm(p_l, u / 2, kamp, kamp2; verbose=verbose)
    #         # Fp_new = -u + 3Fm
    #         # println(Fp, " and ", Fm, " solution: ", u, " vs ", Fp_new, " vs ", Fm_new)
    #     else
    #         u = Ver4.Legrendre(0, wp, angle)
    #         Fp_new = 4 * π / para.me * a_s - u
    #         Fm_new = 0.0
    #     end
    #     if abs(Fp_new - Fp) < eps
    #         break
    #     else
    #         if verbose > 1
    #             println("iter $i: $Fp_new vs $Fp")
    #         end
    #     end
    #     # println(Fm_new, spin_spin)
    #     Fp = mix * Fp_new + (1 - mix) * Fp
    #     Fm = mix * Fm_new + (1 - mix) * Fm
    # end
    # if verbose > 0
    #     println("Self-consistent approach: ")
    #     println("Fs = ", Fp)
    #     println("Fa = ", Fm)
    #     println("u = ", u)
    #     println("4π a_s/m", 4π / para.me * a_s)
    # end
    # return Fp, Fm, u
end

function dKO_dΛ(paras, Λgrid)
    # Fs = [KO(paras[li], lambda, lambda; verbose=0)[1] for lambda in Λgrid]
    dFs = -∂R_∂Λ_exchange(paras, Λgrid; ct=false)[1] ./ (1 .+ ∂R_∂f_exchange(paras, Λgrid; ct=false)[1] ./ paras[1].NF)
    dFm = zero(Λgrid.grid)

    # F0 = -Interp.integrate1D(dFs, Λgrid, [Λgrid[1], Λgrid[end]])
    # Fs = KO(paras[1], Λgrid[1], Λgrid[1]; verbose=0, ct=ct)[1]
    # println(F0, " vs ", Fs[1])
    return dFs, dFm
end

"""
Solve the equation at different energy scale
Find f^+ that with a given u: u = <R^+_exchange(f, Λ)>
"""
function ufp(para, u, kamp=para.kF, kamp2=para.kF; verbose=0, init=-1.0, eps=1e-5)
    println("solving ufp with $u")
    function _u(fs)
        p_l = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs, Fa=0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false)
        wp, wm, angle = Ver4.exchange_interaction(p_l, kamp, kamp2; ct=false, verbose=verbose)
        return Ver4.Legrendre(0, wp, angle) - u
    end

    function _u_f(fs)
        p_l = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs, Fa=0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false)
        wp, wm, angle = Ver4.exchange_interaction_df(p_l, kamp, kamp2; ct=false, verbose=verbose)
        return Ver4.Legrendre(0, wp, angle) / para.NF
    end
    return find_zero((_u, _u_f), init, Roots.Newton(), verbose=(verbose > 0))
end

"""
Solve the equation at different energy scale
Find f^- that with a given u: u = <R^-_exchange(f, Λ)>
"""
function ufm(para, u, kamp=para.kF, kamp2=para.kF; verbose=0, init=-0.0, eps=1e-5)
    println("solving ufm with $u")
    function _u(fa)
        p_l = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=0.0, Fa=fa, order=1, mass2=para.mass2, isDynamic=true, isFock=false)
        wp, wm, angle = Ver4.exchange_interaction(p_l, kamp, kamp2; ct=false, verbose=verbose)
        return Ver4.Legrendre(0, wm, angle) - u
    end

    function _u_f(fa)
        p_l = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fa, Fa=0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false)
        wp, wm, angle = Ver4.exchange_interaction_df(p_l, kamp, kamp2; ct=false, verbose=verbose)
        return Ver4.Legrendre(0, wm, angle) / para.NF
    end
    return find_zero((_u, _u_f), init, Roots.Newton(), verbose=(verbose > 0))
end

function check_charge_instability(para, Λgrid; verbose=0, Fs=0.0, N=100, mix=0.5, eps=1e-5, Fp=zeros(length(Λgrid)))

    as = zeros(length(Λgrid))
    u = zeros(length(Λgrid))
    for li in eachindex(Λgrid)
        lambda = Λgrid[li]
        p_l = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=Fs, Fa=-0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false)
        wp, wm, angle = Ver4.exchange_interaction(p_l, lambda, lambda; ct=false, verbose=verbose)
        u[li] = Ver4.Legrendre(0, wp, angle)
        as[li] = (u[li] + p_l.Fs) * para.me / (4π)
    end

    a, idx = findmax(as)

    println("best: $(a) at the scale $(Λgrid[idx]/para.kF)")
    println(as)
    println(u)
    return idx, as
end

function spin_instability(para, kamp=para.kF, kamp2=para.kF; verbose=0, N=100, mix=0.5, eps=1e-5)

    u = 1.0
    Fs = ufp(para, 1.0, kamp, kamp2; verbose=verbose, init=0.0, eps=eps)
    a_s = (Fs + 1.0) * para.me / (4π)
    return Fs, u, a_s
end

# function check_charge_instability(para, Λgrid; verbose=0, a_s=0.0, N=100, mix=0.5, eps=1e-5, Fp=zeros(length(Λgrid)))
#     Fs = [KO(para, lambda, lambda; verbose=verbose, a_s=a_s, N=N, mix=mix, eps=eps, Fp=Fp[li])[1] for (li, lambda) in enumerate(Λgrid)]
#     limit = [para.NF / UEG.polarKW(q, 0, para) - (para.qTF / q)^2 for q in Λgrid]
#     denorm(q, _Fs) = 1 - UEG.KOinstant(q, para.e0, para.dim, para.mass2, _Fs / para.NF, para.kF) * UEG.polarKW(q, 0, para)
#     d = [denorm(q, Fs[li]) for (li, q) in enumerate(Λgrid)]
#     println(Fs)
#     println(limit)
#     println(d)
#     idx = findall(x -> x < 0.0, d)
#     return idx, Fs
# end

"""
    function ∂Rs_∂Λ_exchange(paras, Λgrid; ct=false)
    
     return N_F ∂<R_{k_1-k_2}>/∂Λ
"""
function ∂R_∂Λ_exchange(paras, Λgrid; ct=false)

    function exchange_s(p, kamp, kamp2)
        wp, wm, angle = Ver4.exchange_interaction(p, kamp, kamp2; ct=ct)
        return Ver4.Legrendre(0, wp, angle)
    end

    function exchange_a(p, kamp, kamp2)
        wp, wm, angle = Ver4.exchange_interaction(p, kamp, kamp2; ct=ct)
        return Ver4.Legrendre(0, wm, angle)
    end

    dFs, dFa = zero(Λgrid), zero(Λgrid)
    for li in eachindex(Λgrid)
        lambda = Λgrid[li]
        p_l = paras[li]
        dFs[li] = central_fdm(5, 1)(λ -> exchange_s(p_l, λ, λ), lambda) #use central finite difference method to calculate the derivative
        dFa[li] = central_fdm(5, 1)(λ -> exchange_a(p_l, λ, λ), lambda) #use central finite difference method to calculate the derivative
    end
    return dFs, dFa
end

"""
    function ∂Rs_∂fs_exchange(paras, Λgrid; ct=false)
    
     return N_F ∂<R_{k_1-k_2}>/∂f
     
"""
function ∂R_∂f_exchange(paras, Λgrid; ct=false)
    dFs, dFa = zero(Λgrid), zero(Λgrid)
    for li in eachindex(Λgrid)
        lambda = Λgrid[li]
        p_l = paras[li]
        wp, wm, angle = Ver4.exchange_interaction_df(p_l, lambda, lambda; ct=ct)
        dFs[li] = Ver4.Legrendre(0, wp, angle)
        dFa[li] = Ver4.Legrendre(0, wm, angle)
    end
    return dFs, dFa
end


function gamma4_treelevel_RG(para, Λgrid; verbose=1, rtol=1e-4, mix=0.9)

    fs = [KO(para, lambda, lambda; verbose=0)[1] for lambda in Λgrid]
    us = deepcopy(fs)

    idx = 1
    while true
        paras = [UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs[li], Fa=0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false) for li in eachindex(Λgrid)]

        ∂R_∂Λ = ∂R_∂Λ_exchange(paras, Λgrid; ct=false)[1]
        ∂R_∂f = ∂R_∂f_exchange(paras, Λgrid; ct=false)[1] / para.NF

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
    # fa = [KO(para, lambda, lambda; verbose=0)[1] for lambda in Λgrid]
    fa = zero(fs)
    u = deepcopy(fs)

    idx = 1
    while true
        paras = [UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs[li], Fa=fa[li], order=1, mass2=para.mass2, isDynamic=true, isFock=false) for li in eachindex(Λgrid)]
        # paras = [UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs[li], Fa=-0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false) for li in eachindex(Λgrid)]

        ∂Rs_∂Λ, ∂Ra_∂Λ = ∂R_∂Λ_exchange(paras, Λgrid; ct=false)
        ∂Rs_∂f, ∂Ra_∂f = ∂R_∂f_exchange(paras, Λgrid; ct=false)
        ∂Rs_∂f ./= para.NF
        ∂Ra_∂f ./= para.NF

        dfs = -(∂Rs_∂Λ .+ 3.0 .* ∂Ra_∂Λ) ./ ∂Rs_∂f ./ 2.0
        dfa = -(∂Rs_∂Λ .- ∂Ra_∂Λ) ./ ∂Ra_∂f ./ 2.0
        # dfs = -(∂Rs_∂Λ .+ 3.0 .* ∂Ra_∂Λ)
        # dfa = -(∂Rs_∂Λ .- ∂Ra_∂Λ)
        du = -(∂Rs_∂Λ .- 3.0 .* ∂Ra_∂Λ)
        fs_new = [-Interp.integrate1D(dfs, Λgrid, [Λgrid[idx], Λgrid[end]]) for idx in eachindex(Λgrid)]
        fa_new = [-Interp.integrate1D(dfa, Λgrid, [Λgrid[idx], Λgrid[end]]) for idx in eachindex(Λgrid)]
        u_new = [-Interp.integrate1D(du, Λgrid, [Λgrid[idx], Λgrid[end]]) for idx in eachindex(Λgrid)]

        max_fs = maximum(abs.((fs .- fs_new)))
        max_fa = maximum(abs.((fa .- fa_new)))
        max_u = maximum(abs.((u .- u_new)))
        if verbose > 0
            println("iteration $(idx) with max_fs = $(max_fs), max_fa = $(max_fa) and max_us = $(max_u)")
        end
        if (max_fs / maximum(abs.(fs)) < rtol) && (max_fa / maximum(abs.(fa)) < rtol) && (max_u / maximum(abs.(u)) < rtol)
            if verbose >= 0
                println("total iteration $(idx)")
            end
            if verbose > 0
                kF_idx = searchsortedfirst(Λgrid, para.kF)
                println("Rs_Λ:", ∂Rs_∂Λ)
                println("Ra_Λ:", ∂Ra_∂Λ)
                println("kF_idx: ", kF_idx, " with ", Λgrid[kF_idx] / para.kF, " kF")
                println("Fs(kF): ", fs[kF_idx])
                println("Fa(kF): ", fa[kF_idx])
                println("U(kF): ", u[kF_idx])
            end
            return fs, fa, u, dfs, dfa, du
        end

        @. fs = fs * (1 - mix) + fs_new * (mix)
        @. fa = fa * (1 - mix) + fa_new * (mix)
        @. u = u * (1 - mix) + u_new * (mix)

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

    rs = [5.0,]
    mass2 = [1e-5,]
    beta = [25.0,]
    order = [1,]
    neval = 1e7

    for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)

        _para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=_order, mass2=_mass2, isDynamic=true)
        Λgrid = CompositeGrid.LogDensedGrid(:gauss, [1.0 * _para.kF, 100 * _para.kF], [_para.kF,], 8, 0.01 * _para.kF, 8)
        fs, fa, u, dfs, dfa, du = gamma4_treelevel_RG_2(_para, Λgrid; verbose=1)
        println(fs)
        println(dfs)

        println(fa)
        println(dfa)

        println(u)
        println(du)

        jldopen("data_f.jld2", "a+") do f
            key = "$(UEG.short(_para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (Λgrid, fs, fa, u, dfs, dfa, du)
        end
    end
end