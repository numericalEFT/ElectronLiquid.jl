"""
In this demo, we analysis the contribution of the exchange KO interaction to the Landau parameters
"""

using Lehmann, GreenFunc, CompositeGrids
using ElectronGas: Polarization
using ElectronGas: SelfEnergy
include("../common/interaction.jl")

function exchange_interaction(Fp, Fm, massratio, Cp=Fp, Cm=Fm)
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]

    Wp = zeros(Float64, length(qs))
    Wm = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi], Wm[qi] = Inter.KO_total(q, 0, para;
            pifunc=Polarization.Polarization0_ZeroTemp_Quasiparticle,
            landaufunc=Inter.landauParameterConst,
            Vinv_Bare=Inter.coulombinv,
            counter_term=Inter.countertermConst,
            Fs=-Fp, Fa=-Fm, Cs=-Cp, Ca=-Cm, massratio=massratio)
        # instantS[qi] = Interaction.coulombinv(q, para)[1]
        # println(q, " -> ", Wp[qi] * NF, ", ", Wm[qi] * NF)
    end
    # Wp *= -NF * z^2 # additional minus sign because the interaction is exchanged
    # Wm *= -NF * z^2
    Wp *= -NF * massratio  # additional minus sign because the interaction is exchanged
    Wm *= -NF * massratio
    return Wp, Wm, θgrid
end

function exchange_interaction_oneloop(Fp, Fm, massratio=1.0, Cp=0.0, Cm=0.0)
    θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    qs = [2 * kF * sin(θ / 2) for θ in θgrid.grid]
    fp = Fs / NF / massratio

    Σ = SelfEnergy.G0W0(para; Euv=100 * para.EF, maxK=8 * para.kF, Nk=16, order=8, minK=1e-8 * para.kF, int_type=:ko_const, Fs=-Fp, Fa=-0.0)
    zz = SelfEnergy.zfactor(Σ)
    z1 = 1 - 1 / zz
    println("z1 = ", z1)

    Wp = zeros(Float64, length(qs))

    for (qi, q) in enumerate(qs)
        # Pi = -NF * massratio * lindhard(q / 2.0 / kF)
        # Wp[qi] = -((KOstatic(q) + fp)^2 - (KOstatic(q))^2) * Pi
        Wp[qi] = 2 * z1 * KOstatic(q)
    end

    Wp *= NF * massratio  # additional minus sign because the interaction is exchanged
    Wm = Wp .* 0.0
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

function projected_exchange_interaction(l, Fp, Fm, massratio, interaction, verbose=1)
    verbose > 0 && println("l=$l:")
    Wse, Wae, θgrid = interaction(Fp, Fm, massratio)
    Wse0 = Legrendre(l, Wse, θgrid)
    Wae0 = Legrendre(l, Wae, θgrid)
    verbose > 1 && println("Wse_l$l=", Wse0)
    verbose > 1 && println("Wae_l$l=", Wae0)

    Ws0, Wa0 = exchange2direct(Wse0, Wae0)
    verbose > 0 && println("Ws_l$l=", Ws0)
    verbose > 0 && println("Wa_l$l=", Wa0)
    return Ws0, Wa0
end

Ws0, Wa0 = projected_exchange_interaction(0, Fs, Fa, massratio, exchange_interaction)
Ws0, Wa0 = projected_exchange_interaction(0, Fs, Fa, massratio, exchange_interaction_oneloop)
exit(0)

if abspath(PROGRAM_FILE) == @__FILE__
    mix = 0.2
    Fp, Fm = Fs, Fa
    for i = 1:100
        # println("iteration: $i")
        global Fp, Fm
        nFs, nFa = projected_exchange_interaction(0, Fp, Fm, massratio, exchange_interaction, 0)
        Fp = Fp * (1 - mix) + nFs * mix
        Fm = Fm * (1 - mix) + nFa * mix
        Fm = 0.0
    end
    println("Self-consistent approach: ")
    println("Fs = ", Fp)
    println("Fa = ", Fm)

    println("Variational approach: ")
    println("parameter Fs     optimized Fs      optimized Fa")
    for Fp in LinRange(0.0, 4.0, 101)
        nFs, nFa = projected_exchange_interaction(0, -Fp, 0.0, massratio, exchange_interaction, 0)
        println("$Fp     $nFs       $nFa")
    end
end