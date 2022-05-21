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

function counterGreen2(Ek, τ, beta, z1, m1, mu1)
    dz = z1
    dm = m1 - z1 / massratio
    dmu = mu1 - z1 * EF

    if τ ≈ 0.0
        g1 = Spectral.kernelFermiT(-1e-8, Ek, β)
    else
        g1 = Spectral.kernelFermiT(τ, Ek, β)
    end
    g2 = green2(Ek, τ, beta)

    return dz * g1 + (dm * Ek - dmu) * g2
end

# double fermi::ThreePhyGreen(double Tau, const momentum &Mom, bool IsFock) {
#   // if tau is exactly zero, set tau=0^-
#   // cout << Tau << endl;

#   double green, Ek;
#   if (Tau == 0.0) {
#     Tau = -1.0e-10;
#   }

#   double s = 1.0;
#   if (Tau < 0.0) {
#     Tau += Para.Beta;
#     s = -s;
#   } else if (Tau >= Para.Beta) {
#     Tau -= Para.Beta;
#     s = -s;
#   }

#   Ek = Mom.squaredNorm(); // bare propagator
#   if (IsFock)
#     Ek += FockSigma(Mom); // Fock diagram dressed propagator
#   else
#     Ek -= Mu_ideal;

#   // double x = Para.Beta * (Ek - Para.Mu) / 2.0;
#   // double y = 2.0 * Tau / Para.Beta - 1.0;
#   // if (x > 100.0)
#   //   green = exp(-x * (y + 1.0));
#   // else if (x < -100.0)
#   //   green = exp(x * (1.0 - y));
#   // else
#   //   green = exp(-x * y) / (2.0 * cosh(x));
#   if (Ek > 0.0) {
#     double Factor = exp(-Para.Beta * Ek);
#     green =
#         exp(-Ek * Tau) / pow(1.0 + Factor, 3.0) *
#         (Tau * Tau / 2.0 -
#          (Para.Beta * Para.Beta / 2.0 + Para.Beta * Tau - Tau * Tau) * Factor +
#          pow(Para.Beta - Tau, 2.0) * Factor * Factor / 2.0);
#   } else {
#     double Factor = exp(Para.Beta * Ek);
#     green =
#         exp(Ek * (Para.Beta - Tau)) / pow(1.0 + Factor, 3.0) *
#         (Tau * Tau / 2.0 * Factor * Factor -
#          (Para.Beta * Para.Beta / 2.0 + Para.Beta * Tau - Tau * Tau) * Factor +
#          pow(Para.Beta - Tau, 2.0) / 2.0);
#   }
#   // cout << green << " vs " << green1 << endl;
#   // if (abs(green - green1) > 1.0e-10) {
#   //   cout << Para.Beta * Ek << endl;
#   //   ABORT("wrong! " << green << ", " << green1 << endl);
#   // }

#   green *= s;

#   if (std::isfinite(green) == false)
#     ABORT("Step:" << Para.Counter << ", Green is too large! Tau=" << Tau
#                   << ", Ek=" << Ek << ", Green=" << green << ", Mom"
#                   << ToString(Mom));
#   // if (std::isnan(green))
#   //   ABORT("Step:" << Para.Counter << ", Green is too large! Tau=" << Tau
#   //                 << ", Ek=" << Ek << ", Green=" << green << ", Mom"
#   //                 << ToString(Mom));
#   return green;
# }