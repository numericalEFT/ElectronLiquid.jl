using ElectronLiquid
using Measurements
using Printf
using JLD2
using PrettyTables

include("z_RG.jl")

const filename = "ver4_PH.jld2"

"""
    function treelevel(para, kgrid, lgrid)
    
    This function calculates the 3-vertex correction to the one-loop 4-vertex.
    
    ∫dΛ [z_1^Λ ∂<R_exchange(f_Λ, Λ)>/∂f + ∂z_1^Λ/∂f <R_exchange(f_Λ, Λ)>]
"""
function vertex3(para, paras, Λgrid, df_dΛ, l)
    # para, kgrid, lgrid, ver4, ver4_df = datatuple
    _Λgrid, z, z1, ∂z1_∂f = z_flow(para, :KO; ngrid=[-1, 0])
    @assert _Λgrid == Λgrid

    Rs_exchange_s = [Ver4.projected_exchange_interaction(l, paras[li], Ver4.exchange_interaction; kamp=Λgrid[li], kamp2=Λgrid[li], ct=true, verbose=0)[1] for li in eachindex(Λgrid)] # counterterm should be on, because the diagrams in DiagMC calcualted with R_q-f, where -f is the counterterm

    Γ3 = z1 .* Rs_exchange_s
    #println(Rs_exchange_s[1])

    ∂Rs_∂fs_exchange_s = -∂Rs_∂fs_exchange(paras, Λgrid; ct=false) ./ 2 # project to the spin-symmetric channel

    ∂Λ3_∂f = z1 .* ∂Rs_∂fs_exchange_s .+ ∂z1_∂f .* Rs_exchange_s
    Γ3_corr = [Interp.integrate1D(∂Λ3_∂f .* df_dΛ, Λgrid, (Λgrid[li], Λgrid[end])) for li in eachindex(Λgrid)]
    println("z1 = ", z1[1], "Rexch = ", Rs_exchange_s[1])

    return Γ3, Γ3_corr
end

# function singularterm(para, kgrid, lgrid)

# end

function process200(para, datatuple)
    # z = CounterTerm.chemicalpotential_renormalization(para.order, sw, dmu)
    key = (2, 0, 0)
    li = 1
    Λgrid, lgrid, _Fs, ver4, ver4_df = datatuple
    # addbare!(datatuple)

    paras = [UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=_Fs[_li], Fa=0.0, order=1, mass2=para.mass2, isDynamic=true, isFock=false) for _li in eachindex(Λgrid)]

    dF_dΛ = dKO_dΛ(paras, Λgrid)[1]
    df_dΛ = dF_dΛ / para.NF

    ######### test if dF_dΛ is correct #########
    F0 = -Interp.integrate1D(dF_dΛ, Λgrid, [Λgrid[1], Λgrid[end]])
    println("grid quality: $F0 vs $(_Fs[1])")
    @assert abs(F0 - _Fs[1]) < 1e-3 "dF_dΛ is not correct! $F0 vs $(_Fs[1])"


    ######## one-loop correction from the Fs flow, only upupdowndown component is needed #######
    ∂F2_∂f = real(ver4_df[(2, 0, 0)][2, li, :]) # index 2 gives the ud component
    F2ud_corr = [Interp.integrate1D(∂F2_∂f .* df_dΛ, Λgrid, [Λgrid[li], Λgrid[end]]) for li in eachindex(Λgrid)]

    ∂F2_∂f = real(ver4_df[(2, 0, 0)][1, li, :]) # index 2 gives the ud component
    F2uu_corr = [Interp.integrate1D(∂F2_∂f .* df_dΛ, Λgrid, [Λgrid[li], Λgrid[end]]) for li in eachindex(Λgrid)]
    println("dF_ud_dΛ corr: ", F2ud_corr[1])
    println("dF_uu_dΛ corr: ", F2uu_corr[1])

    ########### tree-level correction to the one-loop 3-vertex correction ###########
    ver3, ver3_corr = vertex3(para, paras, Λgrid, df_dΛ, lgrid[li])
    ver3 = 2 .* ver3 #right and left 3-vertex correction
    ver3_corr = 2 .* ver3_corr #right and left 3-vertex correction
    println("tree-level correction to the one-loop 3-vertex correction: ", ver3[1], " and ", ver3_corr[1])

    ########## tree-level correction to the one-loop bubble diagram in the exchange channel ##########
    F1bs, F1ba = Ver4.projected_exchange_interaction(0, paras[1], Ver4.exchange_KOcounter; kamp = para.kF, kamp2=para.kF, ct= true, order=1, verbose=0)
    println("tree-level correction to the one-loop bubble diagram in the exchange channel: ", F1bs, ", ", F1ba)

    ############################ one-loop correction  ###########################################
    F2_uu = real(ver4[(2, 0, 0)][1, 1, :])
    F2_ud = real(ver4[(2, 0, 0)][2, 1, :])
    F2s, F2a = Ver4.ud2sa(F2_uu, F2_ud)
    println(F2s, ", ", F2a)

    F1 = _Fs[1]
    F2 = F2s[1] +F1bs - ver3[1] + F2ud_corr[1]

    println(" ------------------------------------------------  ")
    println("   # rs         order 1           order 2    total")
    println("   $(para.rs)     $F1     $F2      $(F1-F2)")
    println(" ------------------------------------------------  ")

    # df = CounterTerm.fromFile("para_wn_1minus0.csv")
    # mu, sw = CounterTerm.getSigma(df, UEG.paraid(para), para.order)
    # dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)

    # println(sw)
    # println(dz)
    # dz =[-0.509, ]
    # dz = [0.0, 0.0]
    # ver4 = CounterTerm.z_renormalization(2, ver4, dz, 2)
    # println(ver4)

    # kF = para.kF
    # printstyled("l=$(lgrid[li]) for $(UEG.short(para))\n", color=:green)
    kgrid = Λgrid[1:1]
    for (p, data) in ver4
        # head = ["k/kF", "uu", "ud", "symmetric", "asymmetric"]
        for (ki, k) in enumerate(kgrid)
            if p[1]>1
                ########## tree-level correction to the one-loop bubble diagram in the exchange channel ##########
                F1bs, F1ba = Ver4.projected_exchange_interaction(lgrid[li], paras[li], Ver4.exchange_KOcounter; kamp = k, kamp2=k, ct= true, order=p[1]-1, verbose=0)
                F1uu, F1ud = Ver4.ud2sa(F1bs, F1ba)
                println("tree-level correction to the one-loop bubble diagram in the exchange channel: ", F1uu, ", ", F1ud)
                println(data[1, li, ki])
                data[1, li, ki] += F1uu - ver3[ki]*2
                data[2, li, ki] += F1ud + F2ud_corr[ki]
            end
        end

        printstyled(@sprintf("%12s    %24s    %24s    %24s    %24s     %24s\n",
                "k/kF", "uu", "ud", "symmetric", "asymmetric", "p: $p"), color=:yellow)
        for (ki, k) in enumerate(kgrid)
            factor = 1.0
            d1, d2 = real(data[1, li, ki]) * factor, real(data[2, li, ki]) * factor
            # s, a = (d1 + d2) / 2.0, (d1 - d2) / 2.0
            s, a = Ver4.ud2sa(d1, d2)
            @printf("%12.6f    %24s    %24s    %24s    %24s\n", k / para.kF, "$d1", "$d2", "$s", "$a")
        end
    end

    printstyled("summed\n", color=:red)
    printstyled(@sprintf("%12s    %24s    %24s    %24s    %24s\n",
            "k/kF", "uu", "ud", "symmetric", "asymmetric"), color=:yellow)
    for (ki, k) in enumerate(kgrid)
        data = ver4[(1, 0)] .+ ver4[(2, 0)]
        d1, d2 = real(data[1, li, ki]), real(data[2, li, ki])
        # s, a = (d1 + d2) / 2.0, (d1 - d2) / 2.0
        s, a = Ver4.ud2sa(d1, d2)
        @printf("%12.6f    %24s    %24s    %24s    %24s\n", k / para.kF, "$d1", "$d2", "$s", "$a")
    end

end

if abspath(PROGRAM_FILE) == @__FILE__

    # @assert length(ARGS) >= 1 "One argument for the data file name is required!"
    # filename = ARGS[1]
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    # rs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
    rs = [4.0, ]
    mass2 = [0.01,]
    _Fs = [-0.0,]
    beta = [25.0,]
    order = [1,]

    f = jldopen(filename, "r")

    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, _Fs, beta, order)
        para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=true)
        kF = para.kF
        for key in keys(f)
            # println(key, " vs ", UEG.paraid(para))
            if key == UEG.short(para)
                process200(para, f[key])
            end
        end
    end

end