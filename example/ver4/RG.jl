using ElectronLiquid
using CompositeGrids
using Plots

rs = 1.0
mass2 = 0.001
beta = 25.0
Nk = 5000

para = UEG.ParaMC(rs=rs, beta=beta, Fs=-0.0, order=1, mass2=mass2, isDynamic=true, isFock=false)
# Λgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, 10 * para.kF], [para.kF,], 8, 0.01 * para.kF, 12)
Λgrid = collect(reverse(LinRange(0.0, 50 * para.kF, Nk)))
dΛgrid = [Λgrid[li+1] - Λgrid[li] for li in 1:length(Λgrid)-1]
println(length(Λgrid))
dΛ = (Λgrid[2] - Λgrid[1]) / 100.0

f = [0.0,]

for (li, lambda) in enumerate(Λgrid[1:end-1])
    f_l = f[end]
    p_l = UEG.ParaMC(rs=rs, beta=beta, Fs=-f_l, order=1, mass2=mass2, isDynamic=true, isFock=false)
    wp, wm, angle = Ver4.exchange_interaction(p_l, p_l.kF, lambda; ct=0.0)
    Fs = Ver4.Legrendre(0, wp, angle) / 2.0

    wp, wm, angle = Ver4.exchange_interaction(p_l, p_l.kF, lambda + dΛ; ct=0.0)
    Fs_dΛ = Ver4.Legrendre(0, wp, angle) / 2.0
    dFs = (Fs_dΛ - Fs) / dΛ
    f_new = f_l + dFs * (Λgrid[li+1] - Λgrid[li])
    push!(f, f_new)
end
kF_idx = searchsortedlast(Λgrid, para.kF, rev=true)
println("kF_idx: ", kF_idx, " with ", Λgrid[kF_idx] / para.kF, " kF")
println("F(kF): ", f[kF_idx-1])

p = plot(Λgrid ./ para.kF, f, label="F", xlabel="Λ", ylabel="F", legend=:bottomright, title="F vs Λ", size=(800, 600), xlims=[0.0, 5.0])
#plot a vertical line at 2kF
plot!(p, [2.0, 2.0], [minimum(f), maximum(f)], label="2k_F", color=:red, linestyle=:dash)
plot!(p, [1.0, 1.0], [minimum(f), maximum(f)], label="k_F", color=:blue, linestyle=:dash)
display(p)
readline()