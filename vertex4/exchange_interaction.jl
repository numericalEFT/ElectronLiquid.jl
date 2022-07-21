include("../common/para_builder.jl")
using .UEG
using CompositeGrids

function exchange_interaction(para::ParaMC)
    kF = para.kF
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]

    Wp = zeros(Float64, length(qs))
    Wm = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = KOstatic(q, para)
        # Wp[qi], Wm[qi] = Inter.KO_total(q, 0, para;
        #     pifunc=Polarization.Polarization0_ZeroTemp_Quasiparticle,
        #     landaufunc=Inter.landauParameterConst,
        #     Vinv_Bare=Inter.coulombinv,
        #     counter_term=Inter.countertermConst,
        #     Fs=-Fp, Fa=-Fm, Cs=-Cp, Ca=-Cm, massratio=massratio)
        # instantS[qi] = Interaction.coulombinv(q, para)[1]
        # println(q, " -> ", Wp[qi] * NF, ", ", Wm[qi] * NF)
    end
    # Wp *= -NF * z^2 # additional minus sign because the interaction is exchanged
    # Wm *= -NF * z^2
    Wp *= -para.NFstar  # additional minus sign because the interaction is exchanged
    Wm *= -para.NFstar
    return Wp, Wm, θgrid
end

function exchange_Coulomb(para::ParaMC)
    kF = para.kF
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]

    Wp = zeros(Float64, length(qs))
    Wm = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = Coulombinstant(q, para)
        # instantS[qi] = Interaction.coulombinv(q, para)[1]
        # println(q, " -> ", Wp[qi] * NF, ", ", Wm[qi] * NF)
    end
    Wp *= -para.NFstar
    return Wp, Wm, θgrid
end

function exchange_KOcounter(para::ParaMC, order, bubble)
    kF = para.kF
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]

    # non-proper and no bubble diagram. (MC for vertex4 doesn't include the buble contribution)
    Wp = UEG.counterKO_W(para; qgrid=qs, ngrid=[0,], order=order, proper=false, bubble=false)[:, 1]
    Wp *= -para.NFstar
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

function projected_exchange_interaction(l, para, interaction, verbose=1)
    verbose > 0 && println(UEG.short(para))
    verbose > 0 && println("l=$l:")
    Wse, Wae, θgrid = interaction(para)
    Wse0 = Legrendre(l, Wse, θgrid)
    Wae0 = Legrendre(l, Wae, θgrid)
    verbose > 1 && println("Wse_l$l=", Wse0)
    verbose > 1 && println("Wae_l$l=", Wae0)

    Ws0, Wa0 = exchange2direct(Wse0, Wae0)
    verbose > 0 && println("Ws_l$l=", Ws0)
    verbose > 0 && println("Wa_l$l=", Wa0)
    verbose > 0 && println()
    return Ws0, Wa0
end

function self_consistent_F0(para::ParaMC, N=100, mix=0.2)
    Fp, Fm = Fs, Fa
    for i = 1:N
        # println("iteration: $i")
        Fp, Fm
        nFs, nFa = projected_exchange_interaction(0, para, exchange_interaction, 0)
        Fp = Fp * (1 - mix) + nFs * mix
        Fm = Fm * (1 - mix) + nFa * mix
        Fm = 0.0
    end
    println("Self-consistent approach: ")
    println("Fs = ", Fp)
    println("Fa = ", Fm)
end

const Order = 1

######################## F1 ##########################
# p = ParaMC(rs=5.0, beta=100.0, Fs=-0.0, order=1, mass2=1e-5)
# projected_exchange_interaction(0, p, exchange_interaction)
# p = ParaMC(rs=5.0, beta=100.0, Fs=-1.0, order=1, mass2=1e-5)
# projected_exchange_interaction(0, p, exchange_interaction)
p = ParaMC(rs=5.0, beta=25.0, Fs=-0.585, order=1, mass2=0.001)
# p = ParaMC(rs=5.0, beta=100.0, Fs=-0.585, order=1, mass2=1e-5)
projected_exchange_interaction(0, p, exchange_interaction)

################# F2 counterterm ##########################
# p = ParaMC(rs=5.0, beta=100.0, Fs=-0.0, order=1, mass2=1e-5)
# projected_exchange_interaction(0, p, p -> exchange_KOcounter(p, 1, false))
# p = ParaMC(rs=5.0, beta=100.0, Fs=-1.0, order=1, mass2=1e-5)
# projected_exchange_interaction(0, p, p -> exchange_KOcounter(p, 1, false))
projected_exchange_interaction(0, p, p -> exchange_KOcounter(p, 1, true))

# println("Variational approach: ")
# println("parameter Fs     optimized Fs      optimized Fa")
# for Fp in LinRange(0.0, 4.0, 101)
#     nFs, nFa = projected_exchange_interaction(0, para, exchange_interaction, 0)
#     println("$Fp     $nFs       $nFa")
# end