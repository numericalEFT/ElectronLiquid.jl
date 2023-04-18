using ElectronLiquid
using CompositeGrids
using Plots
using LaTeXStrings
using FiniteDifferences
pgfplotsx()

include("KO.jl")

rs = 8.0
mass2 = 1e-5
beta = 100.0

para = UEG.ParaMC(rs=rs, beta=beta, Fs=-0.0, order=1, mass2=mass2, isDynamic=true, isFock=false)
Λgrid = CompositeGrid.LogDensedGrid(:gauss, [1.0 * para.kF, 20 * para.kF], [para.kF,], 16, 0.01 * para.kF, 16)
println(length(Λgrid))
println(Λgrid)
dΛ = 0.001 * para.kF

fs = [KO(para, lambda, lambda; verbose=0)[1] for lambda in Λgrid]
us = deepcopy(fs)

dfs = zero(Λgrid.grid)
dus = zero(Λgrid.grid)

function Fs(para, kamp, kamp2)
    wp, wm, angle = Ver4.exchange_interaction(para, kamp, kamp2; ct=false)
    return -Ver4.Legrendre(0, wp, angle)
end

function Fa(para, kamp, kamp2)
    wp, wm, angle = Ver4.exchange_interaction(para, kamp, kamp2; ct=false)
    return -Ver4.Legrendre(0, wm, angle)
end

idx = 1
while true
    println("iteration $(idx)")
    fs_new, us_new = zero(fs), zero(us)
    flag = true
    for (li, lambda) in enumerate(Λgrid)
        p_l = UEG.ParaMC(rs=rs, beta=beta, Fs=fs[li], Fa=0.0, order=1, mass2=mass2, isDynamic=true, isFock=false)

        # _Fs_dΛ = Fs(p_l, lambda + dΛ, lambda + dΛ)
        # _Fs = Fs(p_l, lambda, lambda)
        # dFs = (_Fs_dΛ - _Fs) / dΛ
        dFs = central_fdm(5, 1)(λ -> Fs(p_l, λ, λ), lambda) #use central finite difference method to calculate the derivative

        wp, wm, angle = Ver4.exchange_interaction_df(p_l, lambda, lambda; ct=false)
        Fs_df = Ver4.Legrendre(0, wp, angle) / para.NF

        dfs[li] = -dFs / Fs_df / 2.0
        dus[li] = dFs / 2.0

        fs_new[li] = Interp.integrate1D(dfs, Λgrid, [Λgrid[li], Λgrid[end]])
        us_new[li] = Interp.integrate1D(dus, Λgrid, [Λgrid[li], Λgrid[end]])
        if abs(fs[li] - fs_new[li]) > 1e-4 || abs(us[li] - us_new[li]) > 1e-4
            flag = false
        end
    end

    if flag
        break
    end
    fs .= fs_new
    us .= us_new
    global idx += 1
end
# kF_idx = searchsortedlast(Λgrid, para.kF, rev=true)
kF_idx = 1
println("kF_idx: ", kF_idx, " with ", Λgrid[kF_idx] / para.kF, " kF")
println("Fs(kF): ", fs[kF_idx])
println("Us(kF): ", us[kF_idx])
# println("Fa(kF): ", fa[kF_idx])


############## construct KO interaction ###########
KOgrid = collect(LinRange(1.0 * para.kF, 5.0 * para.kF, 100))
# F_ko = [KO(para, para.kF, lambda; verbose=0)[1] for lambda in KOgrid]
F_ko = [KO(para, lambda, lambda; verbose=0)[1] for lambda in KOgrid]

F_cs = Interp.interp1DGrid(fs, Λgrid, KOgrid)
U_cs = Interp.interp1DGrid(us, Λgrid, KOgrid)

# kF_idx = searchsortedfirst(KOgrid, para.kF)
kF_idx = 1
println("kF_idx: ", kF_idx, " with ", KOgrid[kF_idx] / para.kF, " kF")
println("KO: Fs(kF): ", F_ko[kF_idx])
println(KO(para, para.kF, para.kF))
##################################################

p = plot(titlefontsize=18,
    guidefontsize=18,
    tickfontsize=18,
    legendfontsize=18,
    labelfontsize=18,
    ylabel="Quasiparticle Interaction", xlabel=L"\Lambda/k_F", legend=:topright, size=(800, 600),
    xlims=[0.0, 3.0],
    ylims=[0.0, 1.2],
    linewidth=2, thickness_scaling=1
)
# plot!(p, (Λgrid .- para.kF) ./ para.kF, linewidth=2, fs .* Λgrid, label=L"CS")
# plot!(p, (KOgrid .- para.kF) ./ para.kF, linewidth=2, F_ko .* KOgrid, label=L"KO")

############# plot f ##################
# plot!(p, (Λgrid .- para.kF) ./ para.kF, linewidth=2, fs, label=L"CS")
# plot!(p, (KOgrid .- para.kF) ./ para.kF, linewidth=2, F_ko, label=L"KO")

############ plot u ############
plot!(p, (KOgrid .- para.kF) ./ para.kF, linewidth=2, -F_cs, label=L"-F^{CS}_s")
plot!(p, (KOgrid .- para.kF) ./ para.kF, linewidth=2, U_cs, label=L"U^{CS}_s")
plot!(p, (KOgrid .- para.kF) ./ para.kF, linewidth=2, F_ko, label=L"-F^{KO}_s=U^{KO}_s")
#plot a vertical line at 2kF
# savefig(p, "Us.pdf")
display(p)
readline()