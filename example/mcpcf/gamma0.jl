using ElectronLiquid
using ElectronLiquid.CompositeGrids
using ElectronLiquid.UEG

# rs = 3.0
# beta = 25
# Fs = -0.433
# order = 2
# neval = 1e7

# paramc = UEG.ParaMC(rs=rs, beta=beta, Fs=Fs, order=order, isDynamic=true)
# UEG.MCinitialize!(paramc)

# println(paramc)

function PP_interaction_dynamic(n, para::ParaMC, kamp=para.kF, kamp2=para.kF; ct=true, kawargs...)
    kF = para.kF

    # θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    xgrid = CompositeGrid.LogDensedGrid(:gauss, [-1.0, 1.0], [-1.0, 1.0], 16, 0.001, 16)
    # qs = [2 * kamp * sin(θ / 2) for θ in θgrid.grid]
    # qs = [sqrt(kamp^2 + kamp2^2 - 2 * cos(θ) * kamp * kamp2) for θ in θgrid.grid]
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * x * kamp * kamp2) for x in xgrid]

    Wp = zeros(Float64, length(qs))
    # Wm = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        # Wp[qi] = KOstatic(q, para; ct=ct)
        Wp[qi] = UEG.KO_W(q, n, para)
        # Wp[qi] = UEG.KO_W(q, n, para) * sin(θgrid[qi])
        # Wm[qi] = UEG.KOstatic_spin(q, para; ct=ct)
    end
    Wp *= para.NFstar  # additional minus sign because the interaction is exchanged
    # Wm *= para.NFstar
    # return Wp, Wm, θgrid
    # return Interp.integrate1D(Wp, θgrid)
    return Interp.integrate1D(Wp, xgrid)
end

function gamma1(rslist, nlist)
    results = zeros(Float64, (length(rslist), length(nlist)))
    for irs in 1:length(rslist)
        rs = rslist[irs]
        # paramc = UEG.ParaMC(rs=rs, beta=400, isDynamic=true, Fs=-0.2 * rs)
        paramc = UEG.ParaMC(rs=rs, beta=25, isDynamic=true, Fs=-0.0)
        UEG.MCinitialize!(paramc)
        for iw in 1:length(nlist)
            results[irs, iw] = PP_interaction_dynamic(nlist[iw], paramc)
        end
    end
    return results
end
# println(PP_interaction_dynamic(1, paramc))
rslist = [1.0i for i in 1:3]
nlist = [0, 1]
results = gamma1(rslist, nlist)
println(results)

# using Plots
# p = plot()
# plot!(p, rslist, results[:, 1])
# plot!(p, rslist, results[:, 2])
# plot!(p, rslist, results[:, 3])
# display(p)
# readline()