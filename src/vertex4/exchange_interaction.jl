"""
    function exchange_interaction(para::ParaMC, kamp=para.kF, kamp2=para.kF; ct=true, kwargs...)
    
    Calculate the exchange interaction between two particles with two incoming momentum kamp and kamp2. Averge over the angle between the two momenta.
"""
function exchange_interaction(para::ParaMC, kamp=para.kF, kamp2=para.kF; ct=true,
    θgrid=CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 32),
    kwargs...)
    kF = para.kF
    # println(kamp, ", ", kamp2)
    # qs = [2 * kamp * sin(θ / 2) for θ in θgrid.grid]
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * cos(θ) * kamp * kamp2) for θ in θgrid.grid]

    Wp = zeros(Float64, length(qs))
    Wm = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = KOstatic(q, para; ct=ct)
        Wm[qi] = UEG.KOstatic_spin(q, para; ct=ct)
    end
    # Wp *= -NF * z^2 # additional minus sign because the interaction is exchanged
    # Wm *= -NF * z^2
    Wp *= para.NFstar  # additional minus sign because the interaction is exchanged
    Wm *= para.NFstar
    return Wp, Wm, θgrid
end

function exchange_interaction_df(para::ParaMC, kamp=para.kF, kamp2=para.kF; ct=true, kwargs...)
    kF = para.kF
    # println(kamp, ", ", kamp2)
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    # qs = [2 * kamp * sin(θ / 2) for θ in θgrid.grid]
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * cos(θ) * kamp * kamp2) for θ in θgrid.grid]
    # println(qs)

    Wp = zeros(Float64, length(qs))
    Wm = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = UEG.KOstatic_df(q, para; ct=ct)
        Wm[qi] = UEG.KOstatic_spin_df(q, para; ct=ct)
    end
    # Wp *= -NF * z^2 # additional minus sign because the interaction is exchanged
    # Wm *= -NF * z^2
    Wp *= para.NFstar  # additional minus sign because the interaction is exchanged
    Wm *= para.NFstar
    return Wp, Wm, θgrid
end

"""
    function PP_interaction(para::ParaMC, kamp=para.kF, kamp2=para.kF; ct=true, kawargs...)

 Calculate the interaction averged over the left incoming momentum (0, 0, kamp) and the left outgoing momentum kamp2*(cosθsinϕ, sinθsinϕ, cosϕ).
 The transfer momentum is q = sqrt(kamp^2 + kamp2^2 - 2 * cos(θ) * kamp * kamp2).
"""
function PP_interaction(para::ParaMC, kamp=para.kF, kamp2=para.kF; ct=true, kawargs...)
    kF = para.kF

    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    # qs = [2 * kamp * sin(θ / 2) for θ in θgrid.grid]
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * cos(θ) * kamp * kamp2) for θ in θgrid.grid]

    Wp = zeros(Float64, length(qs))
    Wm = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = KOstatic(q, para; ct=ct)
        Wm[qi] = UEG.KOstatic_spin(q, para; ct=ct)
    end
    Wp *= para.NFstar  # additional minus sign because the interaction is exchanged
    Wm *= para.NFstar
    return Wp, Wm, θgrid
end


function exchange_Coulomb(para::ParaMC, kamp=para.kF, kamp2=para.kF; kwargs...)
    kF = para.kF
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    # qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * cos(θ) * kamp * kamp2) for θ in θgrid.grid]

    Wp = zeros(Float64, length(qs))
    Wm = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = Coulombinstant(q, para)
        # instantS[qi] = Interaction.coulombinv(q, para)[1]
        # println(q, " -> ", Wp[qi] * NF, ", ", Wm[qi] * NF)
    end
    Wp *= para.NFstar
    return Wp, Wm, θgrid
end

function exchange_KOcounter(para::ParaMC, kamp=para.kF, kamp2=para.kF; order, bubble=false, kwargs...)
    kF = para.kF
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    # qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * cos(θ) * kamp * kamp2) for θ in θgrid.grid]

    # non-proper and no bubble diagram. (MC for vertex4 doesn't include the buble contribution)
    Wp = UEG.counterKO_W(para; qgrid=qs, ngrid=[0,], order=order, proper=false, bubble=bubble)[:, 1]
    Wp *= para.NFstar
    Wm = zeros(Float64, length(qs))
    return Wp, Wm, θgrid
end

function Legrendre(l, func, θgrid)
    if l == 0
        return Interp.integrate1D(func .* sin.(θgrid.grid), θgrid) / 2
    elseif l == 1
        return Interp.integrate1D(func .* cos.(θgrid.grid) .* sin.(θgrid.grid), θgrid) / 2
    else
        error("not implemented!")
    end
end

# exchange interaction (Ws + Wa \sigma\sigma)_ex to a direct interaction Ws'+Wa' \sigma\sigma 
# # exchange S/A interaction projected to the spin-symmetric and antisymmetric parts
function exchange2direct(Wse, Wae)
    Ws = (Wse + 3 * Wae) / 2
    Wa = (Wse - Wae) / 2
    return Ws, Wa
end

function projected_exchange_interaction(l, para, interaction; kamp=para.kF, kamp2=para.kF, verbose=1, kwargs...)
    # verbose > 0 && println(UEG.short(para))
    verbose > 0 && println("l=$l:")
    Wse, Wae, θgrid = interaction(para, kamp, kamp2; kwargs...)
    # println(Wse)
    # println(Wae)
    # exit(0)
    Wse0 = Legrendre(l, Wse, θgrid)
    Wae0 = Legrendre(l, Wae, θgrid)
    verbose > 1 && println("Wse_l=$l=", Wse0)
    verbose > 1 && println("Wae_l=$l=", Wae0)

    Ws0, Wa0 = exchange2direct(Wse0, Wae0)
    verbose > 0 && println("Ws_l=$l=", Ws0)
    verbose > 0 && println("Wa_l=$l=", Wa0)
    verbose > 0 && println()
    return Ws0, Wa0
end

function self_consistent_F0(para::ParaMC, kamp=para.kF, kamp2=para.kF, N=100, mix=0.2, verbose=1)
    Fp, Fm = para.Fs, para.Fa
    for i = 1:N
        # println("iteration: $i")
        # Fp, Fm

        para = UEG.reconstruct(para, Fs=-Fp, Fa=-Fm)
        println(para.Fs)
        nFs, nFa = projected_exchange_interaction(0, para, exchange_interaction, verbose=0, kamp=kamp, kamp2=kamp2)
        println("Fp = $nFs, Fm = $nFa")
        Fp = Fp * (1 - mix) + nFs * 2.0 * mix
        # Fm = Fm * (1 - mix) + nFa * mix
        Fm = 0.0
    end
    if verbose > 0
        println("Self-consistent approach: ")
        println("Fs = ", Fp)
        println("Fa = ", Fm)
    end
    return Fp, Fm
end

######################## F1 ##########################
# p = ParaMC(rs=5.0, beta=100.0, Fs=-0.0, order=1, mass2=1e-5)
# projected_exchange_interaction(0, p, exchange_interaction)
# p = ParaMC(rs=5.0, beta=100.0, Fs=-1.0, order=1, mass2=1e-5)
# projected_exchange_interaction(0, p, exchange_interaction)
# p = ParaMC(rs=5.0, beta=25.0, Fs=-0.585, order=1, mass2=0.001)
# p = ParaMC(rs=5.0, beta=100.0, Fs=-0.585, order=1, mass2=1e-5)
# projected_exchange_interaction(0, p, exchange_interaction)

################# F2 counterterm ##########################
# p = ParaMC(rs=5.0, beta=100.0, Fs=-0.0, order=1, mass2=1e-5)
# projected_exchange_interaction(0, p, p -> exchange_KOcounter(p, 1, false))
# p = ParaMC(rs=5.0, beta=100.0, Fs=-1.0, order=1, mass2=1e-5)
# projected_exchange_interaction(0, p, p -> exchange_KOcounter(p, 1, false))
# projected_exchange_interaction(0, p, p -> exchange_KOcounter(p, 1, true))

# println("Variational approach: ")
# println("parameter Fs     optimized Fs      optimized Fa")
# for Fp in LinRange(0.0, 4.0, 101)
#     nFs, nFa = projected_exchange_interaction(0, para, exchange_interaction, 0)
#     println("$Fp     $nFs       $nFa")
# end