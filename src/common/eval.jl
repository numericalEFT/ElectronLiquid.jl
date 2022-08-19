# include("./para_builder.jl")
# using .UEG
module Propagator
using ..FeynmanDiagram
using ..UEG
using ..Lehmann
using LinearAlgebra
using ..ElectronGas

function green2(Ek, τ, beta)
    if τ ≈ 0.0
        τ = -1.0e-10
    end

    s = 1.0
    if τ < 0.0
        τ += beta
        s = -s
    elseif τ >= beta
        τ -= beta
        s = -s
    end

    if Ek > 0.0
        c = exp(-beta * Ek)
        green = exp(-Ek * τ) / (1.0 + c)^2 * (τ - (beta - τ) * c)
    else
        c = exp(beta * Ek)
        green = exp(Ek * (beta - τ)) / (1.0 + c)^2 * (τ * c - (beta - τ))
    end

    green *= s
    #   if (isfinite(green) == false)
    #     ABORT("Step:" << Para.Counter << ", Green is too large! Tau=" << Tau
    #                   << ", Ek=" << Ek << ", Green=" << green << ", Mom"
    #                   << ToString(Mom));
end

function green3(Ek, τ, beta=β)
    if τ ≈ 0.0
        τ = -1.0e-10
    end

    s = 1.0
    if τ < 0.0
        τ += beta
        s = -s
    elseif τ >= beta
        τ -= beta
        s = -s
    end

    if (Ek > 0.0)
        c = exp(-beta * Ek)
        green = exp(-Ek * τ) / (1.0 + c)^3.0 * (τ^2 / 2 - (beta^2 / 2 + beta * τ - τ^2) * c + (beta - τ)^2 * c^2 / 2.0)
    else
        c = exp(beta * Ek)
        green = exp(Ek * (beta - τ)) / (1.0 + c)^3 * (τ^2 * c^2 / 2.0 - (beta^2 / 2.0 + beta * τ - τ^2) * c + (beta - τ)^2 / 2.0)
    end

    green *= s
    return green
end

##################### propagator and interaction evaluation ##############
function DiagTree.eval(id::BareGreenId, K, extT, varT, p::ParaMC)
    kF, β, me, μ, massratio = p.kF, p.β, p.me, p.μ, p.massratio
    τin, τout = varT[id.extT[1]], varT[id.extT[2]]
    k = norm(K)
    if p.isFock
        fock = SelfEnergy.Fock0_ZeroTemp(k, p.basic) - SelfEnergy.Fock0_ZeroTemp(kF, p.basic)
        ϵ = k^2 / (2me * massratio) - μ + fock
    else
        ϵ = k^2 / (2me * massratio) - μ
        # ϵ = kF / me * (k - kF)
    end

    # if k < 0.4 * kF || k > kF * 1.3
    #     return 0.0
    # end

    τ = τout - τin
    order = id.order[1]
    if order == 0
        if τ ≈ 0.0
            return Spectral.kernelFermiT(-1e-8, ϵ, β)
        else
            return Spectral.kernelFermiT(τ, ϵ, β)
        end
    elseif order == 1
        return green2(ϵ, τ, β)
    elseif order == 2
        return green3(ϵ, τ, β)
    elseif order == 3
        return green3(ϵ, τ, β)
    else
        error("not implemented!")
    end
end

# eval(id::InteractionId, K, varT) = e0^2 / ϵ0 / (dot(K, K) + mass2)
function DiagTree.eval(id::BareInteractionId, K, extT, varT, p::ParaMC)
    dim, e0, ϵ0, mass2 = p.dim, p.e0, p.ϵ0, p.mass2
    qd = sqrt(dot(K, K))
    if id.order[2] == 0
        if id.type == Instant
            if interactionTauNum(id.para) == 1
                # return e0^2 / ϵ0 / (dot(K, K) + mass2)
                return Coulombinstant(qd, p)
            elseif interactionTauNum(id.para) == 2
                # println(id.extT)
                return interactionStatic(p, qd, varT[id.extT[1]], varT[id.extT[2]])
            else
                error("not implemented!")
            end
        elseif id.type == Dynamic
            return interactionDynamic(p, qd, varT[id.extT[1]], varT[id.extT[2]])
        else
            error("not implemented!")
        end
    else # counterterm for the interaction
        order = id.order[2]
        if id.type == Instant
            if interactionTauNum(id.para) == 1
                if dim == 3
                    invK = 1.0 / (qd^2 + mass2)
                    return e0^2 / ϵ0 * invK * (mass2 * invK)^order
                    # elseif dim == 2
                    #     invK = 1.0 / sqrt(qd^2 + mass2)
                    #     return e0^2 / ϵ0 * invK * (mass2 * invK)^order
                else
                    error("not implemented!")
                end
            else
                # return counterR(qd, varT[id.extT[1]], varT[id.extT[2]], id.order[2])
                return 0.0 #for dynamical interaction, the counter-interaction is always dynamic!
            end
        elseif id.type == Dynamic
            return counterR(p, qd, varT[id.extT[1]], varT[id.extT[2]], id.order[2])
        else
            error("not implemented!")
        end
    end
end

end