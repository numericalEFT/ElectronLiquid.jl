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

# function counterGreen2(k, τ, beta, mu, me, massratio, z1, m1, mu1)
#     dz = z1
#     # dm = m1 - z1 / massratio
#     dmu = mu1
#     # dmu = mu1
#     Ek = k^2 / (2me * massratio) - mu

#     if τ ≈ 0.0
#         g1 = Spectral.kernelFermiT(-1e-8, Ek, β)
#     else
#         g1 = Spectral.kernelFermiT(τ, Ek, β)
#     end
#     g2 = green2(Ek, τ, beta)

#     # return dz * g1 + (dm * k^2 / (2me) - dmu) * g2
#     # return dz * g1 - dmu * g2  #turn off the mass renormalization
#     return g2  #simply return g*g
#     # return -g2
# end