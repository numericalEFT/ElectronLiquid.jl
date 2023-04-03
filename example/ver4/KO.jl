using ElectronLiquid

"""
    function self_consistent_F0(para::ParaMC, kamp=para.kF, kamp2=para.kF, N=100, mix=0.2, verbose=1)

Equation that solves f = <R_{-k2+k1}> where the averge is taken over the angle between k1 and k2
"""
function KO(para::ParaMC, kamp=para.kF, kamp2=para.kF; N=100, mix=1.0, verbose=1)
    Fp, Fm = -para.Fs, -para.Fa
    for i = 1:N
        p_l = ParaMC(rs=para.rs, beta=para.beta, Fs=-Fp, Fa=-Fm, order=para.order, mass2=para.mass2, isDynamic=false, isFock=false)
        # println(p_l.Fs)
        wp, wm, angle = Ver4.exchange_interaction(p_l, kamp, kamp2; ct=false)
        nFs = Ver4.Legrendre(0, wp, angle)
        nFa = Ver4.Legrendre(0, wm, angle)
        # println("Fp = $nFs, Fold = $(p_l.Fs)")
        # Fp = Fp * (1 - mix) + nFs * mix # 
        # Fm = Fm * (1 - mix) + nFa * mix
        Fm = 0.0
        Fp = nFs
    end
    if verbose > 0
        println("Self-consistent approach: ")
        println("Fs = ", Fp)
        println("Fa = ", Fm)
    end
    return Fp, Fm
end

if abspath(PROGRAM_FILE) == @__FILE__
    p = ParaMC(rs=5.0, beta=25.0, Fs=-0.5, order=1, mass2=1e-5, isDynamic=true, isFock=false)
    KO(p)

    # p = ParaMC(rs=20.0, beta=25.0, Fs=-0.0, order=1, mass2=1e-5, isDynamic=true, isFock=false)
    # self_consistent_F0(p)

    # p = ParaMC(rs=40.0, beta=25.0, Fs=-0.0, order=1, mass2=1e-5, isDynamic=true, isFock=false)
    # self_consistent_F0(p)
end
