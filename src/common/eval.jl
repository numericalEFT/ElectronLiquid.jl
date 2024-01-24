# include("./para_builder.jl")
# using .UEG
module Propagator
using ..FeynmanDiagram
using ..UEG
using ..Lehmann
using LinearAlgebra
using ..ElectronGas

const TAU_CUTOFF = 1e-10

function green_derive(τ, ϵ, β, order)
    if order == 0
        result = green(τ, ϵ, β)
    elseif order == 1
        result = -Spectral.kernelFermiT_dω(τ, ϵ, β)
    elseif order == 2
        result = Spectral.kernelFermiT_dω2(τ, ϵ, β) / 2.0
    elseif order == 3
        result = -Spectral.kernelFermiT_dω3(τ, ϵ, β) / 6.0
    elseif order == 4
        result = Spectral.kernelFermiT_dω4(τ, ϵ, β) / 24.0
    elseif order == 5
        result = -Spectral.kernelFermiT_dω5(τ, ϵ, β) / 120.0
    else
        error("not implemented!")
        # result = Propagator.green(τ, ϵ, β) * 0.0
    end
    return result
end

function interaction_derive(τ1, τ2, K, p::ParaMC, idorder; idtype=Instant, tau_num=1, isLayered=false)
    dim, e0, ϵ0, mass2 = p.dim, p.e0, p.ϵ0, p.mass2
    qd = sqrt(dot(K, K))
    if idorder[2] == 0
        if idtype == Instant
            if tau_num == 1
                # return e0^2 / ϵ0 / (dot(K, K) + mass2)
                if !(isLayered)
                    return Coulombinstant(qd, p)
                else
                    q = sqrt(dot(K, K) + 1e-16)
                    invK = 1.0 / q
                    return e0^2 / 2ϵ0 * invK * tanh(mass2 * q)
                end
            elseif tau_num == 2
                # println(id.extT)
                if idorder[3] == 0
                    return interactionStatic(p, qd, τ1, τ2)
                else # return dR/df for the RG purpose. The static part is zero
                    return 0.0
                end
            else
                error("not implemented!")
            end
        elseif idtype == Dynamic
            if idorder[3] == 0
                return interactionDynamic(p, qd, τ1, τ2)
            else # return dR/df for the RG purpose.
                return UEG.interactionDynamic_df(p, qd, τ1, τ2)
            end
        else
            error("not implemented!")
        end
    else # counterterm for the interaction
        order = idorder[2]
        if idtype == Instant
            if tau_num == 1
                if dim == 3
                    invK = 1.0 / (qd^2 + mass2)
                    return e0^2 / ϵ0 * invK * (mass2 * invK)^order
                elseif dim == 2
                    if !(isLayered)
                        invK = 1.0 / (qd + mass2)
                        return e0^2 / 2ϵ0 * invK * (mass2 * invK)^order
                    else
                        return 0.0
                    end
                else
                    error("not implemented!")
                end
            else
                # return counterR(qd, τ1, τ2, idorder[2])
                return 0.0 #for dynamical interaction, the counter-interaction is always dynamic!
            end
        elseif idtype == Dynamic
            if idorder[3] == 0
                return counterR(p, qd, τ1, τ2, idorder[2])
            else
                return counterR_df(p, qd, τ1, τ2, idorder[2])
            end
        else
            error("not implemented!")
        end
    end
end

function green(τ::T, ω::T, β::T) where {T}
    #generate green function of fermion
    if τ ≈ T(0.0)
        τ = -TAU_CUTOFF
    end
    if τ > T(0.0)
        return ω > T(0.0) ?
               exp(-ω * τ) / (1 + exp(-ω * β)) :
               exp(ω * (β - τ)) / (1 + exp(ω * β))
    else
        return ω > T(0.0) ?
               -exp(-ω * (τ + β)) / (1 + exp(-ω * β)) :
               -exp(-ω * τ) / (1 + exp(ω * β))
    end
end

function green2(Ek, τ, beta)
    if τ ≈ 0.0
        τ = -TAU_CUTOFF
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
        τ = -TAU_CUTOFF
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
        green = exp(-Ek * τ) / (1.0 + c)^3.0 * (τ^2 - (beta^2 + 2 * beta * τ - 2 * τ^2) * c + (beta - τ)^2 * c^2)
    else
        c = exp(beta * Ek)
        green = exp(Ek * (beta - τ)) / (1.0 + c)^3 * (τ^2 * c^2 - (beta^2 + 2 * beta * τ - 2 * τ^2) * c + (beta - τ)^2)
    end

    green *= s
    return green
end

function dispersion(k, p::ParaMC)
    kF, me, μ, massratio = p.kF, p.me, p.μ, p.massratio
    if p.isFock
        fock = SelfEnergy.Fock0_ZeroTemp(k, p.basic) - SelfEnergy.Fock0_ZeroTemp(kF, p.basic)
        ϵ = k^2 / (2me * massratio) - μ + fock
    else
        ϵ = k^2 / (2me * massratio) - μ
        # ϵ = kF / me * (k - kF)
    end
    return ϵ
end

##################### propagator and interaction evaluation ##############
function DiagTree.eval(id::BareGreenId, K, extT, varT, p::ParaMC)
    β = p.β
    τin, τout = varT[id.extT[1]], varT[id.extT[2]]
    k = norm(K)
    ϵ = dispersion(k, p)
    # if k < 0.4 * kF || k > kF * 1.3
    #     return 0.0
    # end

    # Evaluate μ counterterms; note that we have:
    # ∂^n_μ g(ϵₖ - μ, τ) = (-1)^n ∂^n_ω g(ω, τ)
    τ = τout - τin
    order = id.order[1]
    return green_derive(τ, ϵ, β, order) * factorial(order)
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
                if id.order[3] == 0
                    return interactionStatic(p, qd, varT[id.extT[1]], varT[id.extT[2]])
                else # return dR/df for the RG purpose. The static part is zero
                    return 0.0
                end
            else
                error("not implemented!")
            end
        elseif id.type == Dynamic
            if id.order[3] == 0
                return interactionDynamic(p, qd, varT[id.extT[1]], varT[id.extT[2]])
            else # return dR/df for the RG purpose.
                return UEG.interactionDynamic_df(p, qd, varT[id.extT[1]], varT[id.extT[2]])
            end
        else
            error("not implemented!")
        end
    else # counterterm for the interaction
        order = id.order[2]
        if id.type == Instant
            if interactionTauNum(id.para) == 1
                if dim == 3
                    invK = 1.0 / (qd^2 + mass2)
                    return factorial(order) * e0^2 / ϵ0 * invK * (mass2 * invK)^order
                elseif dim == 2
                    invK = 1.0 / (qd + mass2)
                    return factorial(order) * e0^2 / 2ϵ0 * invK * (mass2 * invK)^order
                    # elseif dim == 2
                    #     invK = 1.0 / sqrt(qd^2 + mass2)
                    #     return factorial(order) * e0^2 / ϵ0 * invK * (mass2 * invK)^order
                else
                    error("not implemented!")
                end
            else
                # return counterR(qd, varT[id.extT[1]], varT[id.extT[2]], id.order[2])
                return 0.0 #for dynamical interaction, the counter-interaction is always dynamic!
            end
        elseif id.type == Dynamic
            if id.order[3] == 0
                return factorial(order) * counterR(p, qd, varT[id.extT[1]], varT[id.extT[2]], id.order[2])
            else
                return factorial(order) * counterR_df(p, qd, varT[id.extT[1]], varT[id.extT[2]], id.order[2])
            end
        else
            error("not implemented!")
        end
    end
end

end
