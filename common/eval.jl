##################### propagator and interaction evaluation ##############
function eval(id::BareGreenId, K, extT, varT)
    τin, τout = varT[id.extT[1]], varT[id.extT[2]]
    k = norm(K)
    if isFock
        fock = SelfEnergy.Fock0_ZeroTemp(k, para) - SelfEnergy.Fock0_ZeroTemp(kF, para)
        ϵ = k^2 / (2me * massratio) - μ + fock
    else
        ϵ = k^2 / (2me * massratio) - μ
    end
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
function eval(id::BareInteractionId, K, extT, varT)
    qd = sqrt(dot(K, K))
    if id.order[2] == 0
        if id.type == Instant
            if id.para.interactionTauNum == 1
                # return e0^2 / ϵ0 / (dot(K, K) + mass2)
                return Coulombinstant(qd)
            elseif id.para.interactionTauNum == 2
                return interactionStatic(qd, varT[id.extT[1]], varT[id.extT[2]])
            else
                error("not implemented!")
            end
        elseif id.type == Dynamic
            return interactionDynamic(qd, varT[id.extT[1]], varT[id.extT[2]])
        else
            error("not implemented!")
        end
    else # counterterm for the interaction
        order = id.order[2]
        if id.type == Instant
            if id.para.interactionTauNum == 1
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
            return counterR(qd, varT[id.extT[1]], varT[id.extT[2]], id.order[2])
        end
    end
end