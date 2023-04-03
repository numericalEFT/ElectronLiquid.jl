using ElectronLiquid
using CompositeGrids
using Plots
using LaTeXStrings
pgfplotsx()

include("KO.jl")

rs = 4.0
mass2 = 1e-5
beta = 100.0
Nk = 40000

para = UEG.ParaMC(rs=rs, beta=beta, Fs=-0.0, order=1, mass2=mass2, isDynamic=true, isFock=false)
# Λgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, 10 * para.kF], [para.kF,], 8, 0.01 * para.kF, 12)
Λgrid = collect(reverse(LinRange(para.kF, 20 * para.kF, Nk)))
dΛgrid = [Λgrid[li+1] - Λgrid[li] for li in 1:length(Λgrid)-1]
println(length(Λgrid))
# dΛ = (Λgrid[2] - Λgrid[1]) / 100.0
dΛ = (Λgrid[2] - Λgrid[1])

fs = [0.0,]
fa = [0.0,]

for (li, lambda) in enumerate(Λgrid[1:end-1])
    fs_l, fa_l = fs[end], fa[end]
    # p_l = UEG.ParaMC(rs=rs, beta=beta, Fs=-fs_l, Fa=-fa_l, order=1, mass2=mass2, isDynamic=true, isFock=false)
    p_l = UEG.ParaMC(rs=rs, beta=beta, Fs=-fs_l, Fa=0.0, order=1, mass2=mass2, isDynamic=true, isFock=false)
    # ws, wa = Ver4.projected_exchange_interaction(0, p_l, Ver4.exchange_interaction; kamp=p_l.kF, kamp2=lambda, ct=false, verbose=0)
    # ws_d, wa_d = Ver4.projected_exchange_interaction(0, p_l, Ver4.exchange_interaction; kamp=p_l.kF, kamp2=lambda + dΛ, ct=false, verbose=0)

    # wp, wm, angle = Ver4.exchange_interaction(p_l, p_l.kF, lambda; ct=false)
    wp, wm, angle = Ver4.exchange_interaction(p_l, lambda, lambda; ct=false)
    Fs = Ver4.Legrendre(0, wp, angle)
    Fa = Ver4.Legrendre(0, wm, angle)

    # wp, wm, angle = Ver4.exchange_interaction(p_l, p_l.kF, lambda + dΛ; ct=false)
    wp, wm, angle = Ver4.exchange_interaction(p_l, lambda + dΛ, lambda + dΛ; ct=false)
    Fs_dΛ = Ver4.Legrendre(0, wp, angle)
    Fa_dΛ = Ver4.Legrendre(0, wm, angle)
    dFs = (Fs_dΛ - Fs) / dΛ
    dFa = (Fa_dΛ - Fa) / dΛ

    # wp, wm, angle = Ver4.exchange_interaction_df(p_l, p_l.kF, lambda; ct=false)
    wp, wm, angle = Ver4.exchange_interaction_df(p_l, lambda, lambda; ct=false)
    Fs_df = Ver4.Legrendre(0, wp, angle) / para.NF
    Fa_df = Ver4.Legrendre(0, wm, angle) / para.NF

    # println(lambda, "->", Fs_df)
    # dFs = (ws_d - ws) / dΛ
    # dFa = (wa_d - wa) / dΛ
    fs_new = fs_l + dFs / Fs_df / 2.0 * (Λgrid[li+1] - Λgrid[li])
    # fs_new = fs_l + dFs * (Λgrid[li+1] - Λgrid[li])
    fa_new = fa_l + dFa / Fa_df / 2.0 * (Λgrid[li+1] - Λgrid[li])
    push!(fs, fs_new)
    push!(fa, fa_new)
end
# kF_idx = searchsortedlast(Λgrid, para.kF, rev=true)
kF_idx = length(Λgrid)
println("kF_idx: ", kF_idx, " with ", Λgrid[kF_idx] / para.kF, " kF")
println("Fs(kF): ", fs[kF_idx])
println("Fa(kF): ", fa[kF_idx])


############## construct KO interaction ###########
KOgrid = collect(LinRange(1.0 * para.kF, 5.0 * para.kF, 100))
# F_ko = [KO(para, para.kF, lambda; verbose=0)[1] for lambda in KOgrid]
F_ko = [KO(para, lambda, lambda; verbose=0)[1] for lambda in KOgrid]

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
    ylabel=L"F_s", xlabel=L"\Lambda/k_F", legend=:topright, size=(800, 600), xlims=[0.0, 3.0], ylims=[0.0, 1.2], linewidth=2, thickness_scaling=1
)
plot!(p, (Λgrid .- para.kF) ./ para.kF, linewidth=2, fs, label=L"CS")
plot!(p, (KOgrid .- para.kF) ./ para.kF, linewidth=2, F_ko, label=L"KO")
#plot a vertical line at 2kF
# savefig(p, "Fs.pdf")
# display(p)
# readline()