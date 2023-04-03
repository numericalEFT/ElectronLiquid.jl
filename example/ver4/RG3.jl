using ElectronLiquid
using CompositeGrids
using Plots

rs = 5.0
mass2 = 0.001
beta = 25.0
Nk = 5000

para = UEG.ParaMC(rs=rs, beta=beta, Fs=-0.0, order=1, mass2=mass2, isDynamic=true, isFock=false)
# Λgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, 10 * para.kF], [para.kF,], 8, 0.01 * para.kF, 12)
Λgrid = collect(reverse(LinRange(0.0, 50 * para.kF, Nk)))
dΛgrid = [Λgrid[li+1] - Λgrid[li] for li in 1:length(Λgrid)-1]
println(length(Λgrid))
dΛ = (Λgrid[2] - Λgrid[1]) / 100.0

fs = [0.0,]
fa = [0.0,]

for (li, lambda) in enumerate(Λgrid[1:end-1])
    fs_l, fa_l = fs[end], fa[end]
    # p_l = UEG.ParaMC(rs=rs, beta=beta, Fs=-fs_l, Fa=-fa_l, order=1, mass2=mass2, isDynamic=true, isFock=false)
    p_l = UEG.ParaMC(rs=rs, beta=beta, Fs=0.0, Fa=-fa_l, order=1, mass2=mass2, isDynamic=true, isFock=false)
    # ws, wa = Ver4.projected_exchange_interaction(0, p_l, Ver4.exchange_interaction; kamp=p_l.kF, kamp2=lambda, ct=false, verbose=0)
    # ws_d, wa_d = Ver4.projected_exchange_interaction(0, p_l, Ver4.exchange_interaction; kamp=p_l.kF, kamp2=lambda + dΛ, ct=false, verbose=0)

    wp, wm, angle = Ver4.exchange_interaction(p_l, p_l.kF, lambda; ct=false)
    Fs = Ver4.Legrendre(0, wp, angle)
    Fa = Ver4.Legrendre(0, wm, angle)

    wp, wm, angle = Ver4.exchange_interaction(p_l, p_l.kF, lambda + dΛ; ct=false)
    Fs_dΛ = Ver4.Legrendre(0, wp, angle)
    Fa_dΛ = Ver4.Legrendre(0, wm, angle)
    dFs = (Fs_dΛ - Fs) / dΛ
    dFa = (Fa_dΛ - Fa) / dΛ

    # dFs = (ws_d - ws) / dΛ
    # dFa = (wa_d - wa) / dΛ
    fs_new = fs_l + dFs * (Λgrid[li+1] - Λgrid[li])
    fa_new = fa_l + dFa * (Λgrid[li+1] - Λgrid[li])
    push!(fs, fs_new)
    push!(fa, fa_new)
end
kF_idx = searchsortedlast(Λgrid, para.kF, rev=true)
println("kF_idx: ", kF_idx, " with ", Λgrid[kF_idx] / para.kF, " kF")
println("Fs(kF): ", fs[kF_idx-1])
println("Fa(kF): ", fa[kF_idx-1])

p = plot(Λgrid ./ para.kF, fs, label="Fs", xlabel="Λ", ylabel="F", legend=:bottomright, title="Fs vs Λ", size=(800, 600), xlims=[0.0, 5.0])
plot!(p, Λgrid ./ para.kF, fa, label="Fa", xlabel="Λ", ylabel="F", legend=:bottomright, title="Fa vs Λ", size=(800, 600), xlims=[0.0, 5.0])
#plot a vertical line at 2kF
plot!(p, [2.0, 2.0], [minimum(fs), maximum(fs)], label="2k_F", color=:red, linestyle=:dash)
plot!(p, [1.0, 1.0], [minimum(fs), maximum(fs)], label="k_F", color=:blue, linestyle=:dash)
display(p)
readline()