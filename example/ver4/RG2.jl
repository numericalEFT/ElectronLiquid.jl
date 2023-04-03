using ElectronLiquid
using CompositeGrids
using Roots
using Plots

rs = 5.0
mass2 = 0.001
beta = 25.0
Nk = 100

para = UEG.ParaMC(rs=rs, beta=beta, Fs=-0.0, order=1, mass2=mass2, isDynamic=true, isFock=false)
# Λgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, 10 * para.kF], [para.kF,], 8, 0.01 * para.kF, 12)
Λgrid = collect(reverse(LinRange(0 * para.kF, 3 * para.kF, Nk)))
dΛgrid = [Λgrid[li+1] - Λgrid[li] for li in 1:length(Λgrid)-1]
println(length(Λgrid))
dΛ = (Λgrid[2] - Λgrid[1]) / 100.0

fs = []
fa = []

function dRex(para, Λ, fs, fa, kamp=para.kF, dΛ=kamp / 10000)
    p = UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=-fs, Fa=-fa, order=1, mass2=para.mass2, isDynamic=para.isDynamic, isFock=para.isFock)
    wp, wm, angle = Ver4.exchange_interaction(p, kamp, Λ; ct=false)
    Fs = Ver4.Legrendre(0, wp, angle)
    Fa = Ver4.Legrendre(0, wm, angle)

    wp, wm, angle = Ver4.exchange_interaction(p, kamp, Λ + dΛ; ct=false)
    Fs2 = Ver4.Legrendre(0, wp, angle)
    Fa2 = Ver4.Legrendre(0, wm, angle)
    dFs = (Fs2 - Fs) / dΛ
    dFa = (Fa2 - Fa) / dΛ
    return dFs, dFa
end

for (li, lambda) in enumerate(Λgrid[1:end])
    println(li)
    # fs_l, fa_l = fs[end], fa[end]
    # fs_new = find_zero(x -> dRex(para, lambda, x, 0.0)[1], (-2.0, 2.0), xatol=1e-2, Bisection(), verbose=true)
    # fs_new = find_zero(x -> dRex(para, lambda, x, 0.0)[1], 1.0, atol=1e-2, verbose=true)
    _fs = [dRex(para, lambda, x, 0.0)[1] for x in LinRange(0.0, 2.0, 100)]
    # find the index where _fs[i] and _fs[i+1] have different signs
    fs_new = 0.0
    for i in 1:length(_fs)-1
        if _fs[i] * _fs[i+1] < 0
            fs_new = (i - 1) / 100.0 - 2.0
            # break
        end
    end

    println(fs_new)
    # fs_new = 
    fa_new = 0.0

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