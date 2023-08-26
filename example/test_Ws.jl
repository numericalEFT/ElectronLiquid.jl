using Test
using ElectronLiquid
using FeynmanDiagram
using FiniteDifferences
using Lehmann
using Measurements
using ElectronLiquid
using ElectronLiquid.CompositeGrids
using ElectronLiquid.UEG

using ElectronLiquid.ElectronGas
using ElectronLiquid.ElectronGas.GreenFunc

function DCkernel(param, nlist)
    kF = param.kF
    Euv, rtol = 1000 * param.EF, 1e-10
    Nk, minK, order, maxK = 8, 1e-7kF, 8, 10kF
    int_type = :rpa

    #--- prepare kernel ---
    W = LegendreInteraction.DCKernel0(param;
        Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order,
        int_type=int_type)
    kgrid = W.kgrid
    qgrids = W.qgrids
    ik = locate(kgrid, kF)
    iq = locate(qgrids[ik], kF)
    Wfreq = W.kernel[ik, iq, :] .+ W.kernel_bare[ik, iq]
    dlr = W.dlrGrid
    Wdlr = Lehmann.matfreq2dlr(dlr, Wfreq)
    Wn = Lehmann.dlr2matfreq(dlr, Wdlr, nlist)
    return Wn ./ kF^2
end

@inline function Vinstant(q, para::ParaMC)
    if abs(q) < 1e-16
        q = 1e-16
    end
    return 4π * para.e0^2 / q^2
end

function KO_Ws_reg(n::Integer, para::ParaMC; Pi=UEG.polarKW(1e-16, n, para))
    q = 1e-16
    V = Vinstant(q, para)
    return 1 / (1 - Pi * V)
end

function PPs(n, para::ParaMC, kamp=para.kF, kamp2=para.kF; kawargs...)
    kF = para.kF
    xgrid = CompositeGrid.LogDensedGrid(:gauss, [-1.0, 1.0], [-1.0, 1.0], 32, 1e-12, 8)
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * x * kamp * kamp2) for x in xgrid]

    Wp = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = KO_Ws_reg(n, para) * UEG.Coulombinstant(q, para)
    end
    Wp *= para.NFstar  # additional minus sign because the interaction is exchanged
    return Interp.integrate1D(Wp, xgrid)
end

function PP_interaction_dynamic(n, para::ParaMC, kamp=para.kF, kamp2=para.kF; kawargs...)
    kF = para.kF
    xgrid = CompositeGrid.LogDensedGrid(:gauss, [-1.0, 1.0], [-1.0, 1.0], 32, 1e-12, 8)
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * x * kamp * kamp2) for x in xgrid]

    Wp = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = UEG.KO_W(q, n, para)
    end
    Wp *= para.NFstar  # additional minus sign because the interaction is exchanged
    return Interp.integrate1D(Wp, xgrid)
end

# @testset "PP" begin
seed = 1234
p = (1, 0, 0)
rs = 3.0
Fs = 0.56
# rs = 1.0
# Fs = 0.209
beta = 400
mass2 = 1e-12
neval = 1e6
para = ElectronLiquid.ParaMC(rs=rs, beta=beta, Fs=Fs, order=1, mass2=mass2, isDynamic=true)
UEG.MCinitialize!(para)

# nlist = [0, 1, 2, 3, 10, 100, 2000, 40000]
nlist = [100, 200, 400, 800, 2000, 40000]
W1s = Float64[]
W2s = Float64[]
for i in nlist
    W1 = PP_interaction_dynamic(i, para) / para.NFstar * para.EF
    W2 = PPs(i, para) / para.NFstar * para.EF
    # println(W1)
    push!(W1s, W1)
    push!(W2s, W2)
end
println([π * (2i) / para.β for i in nlist])
println(W1s)
# println(DCkernel(para.basic, nlist))
println(W2s)

# nmax = floor(Int, 0.1para.EF / 2 / π * para.β + 1.5)

# nlist = [i for i in 0:nmax]
# wlist = [π * (2i) / para.β / para.EF for i in nlist]
# # Wr = [(PP_interaction_dynamic(i, para) - PPs(i, para)) / para.NFstar for i in nlist]
# Wr = [(PP_interaction_dynamic(i, para) - PPs(i, para)) for i in nlist]

# nfull = [i for i in 0:500]
# wfull = [π * (2i) / para.β / para.EF for i in nfull]
# Wrfull = [(PP_interaction_dynamic(i, para) - PPs(i, para)) for i in nfull]
# Wfull = [(PP_interaction_dynamic(i, para)) for i in nfull]
# Γr = zeros(Float64, (length(nlist), length(nlist)))
# for i in 1:length(nlist)
#     for j in 1:length(nlist)
#         np, nm = (i + j + 1), abs(i - j)
#         Γr[i, j] += (PP_interaction_dynamic(np, para) - PPs(np, para))
#         Γr[i, j] += (PP_interaction_dynamic(nm, para) - PPs(nm, para))
#     end
# end

# using Plots, LaTeXStrings
# # plt = plot(legend=:none, xlimit=(0.0, 0.15), ylimit=(0.0, Wr[end]))
# plt = plot(legend=:none, xlimit=(0.0, 0.11), ylimit=(0.0, 0.11), zlimit=(0.0, 1.1))
# # plot!(plt, wlist, Wr)
# # # plot!(wfull, Wrfull, inset=(1, bbox(0.05, 0.05, 0.5, 0.25, :bottom, :right)))
# surface!(plt, wlist, wlist, Γr)
# wireframe!(plt, wlist, wlist, Γr)
# xlabel!(plt, L"$\omega_1/E_F$")
# ylabel!(plt, L"$\omega_2/E_F$")
# zlabel!(plt, L"$\Gamma_r^e(\omega_1,\omega_2)N_F$")
# # display(plt)
# # readline()
# savefig(plt, "run/gamma03d.pdf")

# plt = plot(xlimit=(0.0, 7.0))
# plot!(wfull, Wrfull, label=L"$\Gamma_r^e$")
# plot!(wfull, Wfull, label=L"$\Gamma^e$")
# xlabel!(plt, L"$\omega_1/E_F$")
# ylabel!(plt, L"$\Gamma(\omega)N_F$")
# savefig(plt, "run/gamma0full.pdf")
# # end