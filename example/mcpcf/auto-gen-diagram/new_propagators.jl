# provide and test containers of propagators for MC calculations.

module Propagators

using FeynmanDiagram

using ElectronLiquid
using ElectronLiquid.ElectronGas
using ElectronLiquid.ElectronGas.Interaction: RPAwrapped
using ElectronLiquid.ElectronGas.Parameter
using ElectronLiquid.ElectronGas.GreenFunc
using ElectronLiquid.ElectronGas.CompositeGrids
using JLD2

export rpa, interaction, G0, initR, response, coulomb, fake
export calcF!, responsef
export diagram_gen, diagram_eval

# Ri excludes the source term

const rtol = 1e-10 # rtol of DLR 
const nlog_factor = 3.0 # factor controlling how many point per order of magnitude
const βEUV = 1e6 # β*Euv for Rt

# wrapper of functions and parameters
struct Funcs{P,II,IT,RI,RT,FT}
    param::P
    inti::II
    intt::IT
    # Ri::RI
    Ris::RI # store accumulated Ri for self consistent
    Rt::RT # container collecting Rt result, use DLR τ grid 
    # Rtd::RTD # container for mc sampling, use dense τ grid
    fake::II
    Ft::FT # store F with dense grid
end

function f_ake(ii, it)
    # multiplier for better cancellation between V and W
    # V(q)=\int_0^β V(q)fake(q, t) dt
    idlr = to_dlr(it)
    iw = to_imfreq(idlr)
    fi = deepcopy(ii)
    fi.data[1, 1, :] .= iw[1, 1, :] .+ 1.0
    return fi
end

# Funcs(p, ii, it, ri, rt, rtd) = Funcs(p, ii, it, ri, deepcopy(ri), rt, rtd, f_ake(ii, it), rt)
# Funcs(p, ii, it, ri, rtd) = Funcs(p, ii, it, ri, deepcopy(ri), rtd, rtd, f_ake(ii, it), rt)
# Funcs(p, ii, it, ri, ris, rt, rtd) = Funcs(p, ii, it, ri, ris, rt, rtd, f_ake(ii, it), rt)
Funcs(p, ii, it, ri, rt, ft) = Funcs(p, ii, it, ri, rt, f_ake(ii, it), ft)

interaction(k, funcs::Funcs) = interaction(k, funcs.inti)
interaction(t, k, funcs::Funcs) = interaction(t, k, funcs.intt)
response(k, funcs::Funcs; norm=1) = response(k, funcs.Ri; norm=norm)
# response(t, k, funcs::Funcs; norm=1) = response(t, k, funcs.Rt; norm=norm)
response(t, k, funcs::Funcs; norm=1) = response(t, k, funcs.Rtd; norm=norm)
responsef(t, k, funcs::Funcs; norm=1) = response(t, k, funcs.Ft; norm=norm)
G0(t, k, funcs::Funcs) = G0(t, k, funcs.param)
coulomb(q, funcs::Funcs) = coulomb(q, funcs.param)
fake(k, funcs::Funcs) = interaction(k, funcs.fake)

# shift tau to [0, β)
function tau_fermi(t, β)
    if 0 <= t < β
        return t, 1.0
    elseif β <= t < 2β
        return t - β, -1.0
    elseif -β <= t < 0
        return t + β, -1.0
    else
        error("τ=$t out of range!")
    end
end

function tau_bose(t, β)
    if 0 <= t < β
        return t, 1.0
    elseif β <= t < 2β
        return t - β, 1.0
    elseif -β <= t < 0
        return t + β, 1.0
    else
        error("τ=$t out of range!")
    end
end

function rpa(param;
    mint=0.001,
    minK=0.001 * sqrt(param.T * param.me), maxK=10.0 * param.kF,
    order=6)

    Nlog = floor(Int, 2.0 * log10(param.β / mint))
    tgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, param.β], [0.0, param.β], Nlog, mint, order)

    Nk = floor(Int, 2.0 * log10(maxK / minK))
    kgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)

    rpad, rpai = Propagators.RPAwrapped(100 * param.EF, 1e-10, kgrid, param)
    rpadlr = Propagators.to_dlr(rpad)

    rpat = Propagators.dlr_to_imtime(rpadlr, tgrid)
    return rpai, rpat
end

function phonon(param;
    mint=0.001,
    minK=0.001 * sqrt(param.T * param.me), maxK=10.0 * param.kF,
    order=6)

    Nlog = floor(Int, 2.0 * log10(param.β / mint))
    tgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, param.β], [0.0, param.β], Nlog, mint, order)

    Nk = floor(Int, 2.0 * log10(maxK / minK))
    kgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)

    Euv = 100 * param.EF
    rtol = 1e-10
    β = param.β

    wn_mesh = GreenFunc.ImFreq(β, BOSON; Euv=Euv, rtol=rtol, symmetry=:ph)
    green_dyn = GreenFunc.MeshArray(1:2, wn_mesh, kgrid; dtype=ComplexF64)

    for (ki, k) in enumerate(kgrid)
        for (ni, n) in enumerate(wn_mesh.grid)
            green_dyn[1, ni, ki], green_dyn[2, ni, ki] = Interaction.phonon(k, n, param)
        end
    end
    phonondlr = Propagators.to_dlr(green_dyn)
    phonont = Propagators.dlr_to_imtime(phonondlr, tgrid)

    return phonont
end

function coulomb(q, param)
    return ElectronGas.Interaction.coulomb(q, param)
end

function interaction(k, prop; i=1)
    if k ≈ 0
        k = 1e-6
    end
    # return Interp.interp1D(view(prop.data, i, 1, :), prop.mesh[3], k)
    return Interp.linear1D(view(prop.data, i, 1, :), prop.mesh[3], k)
end

function interaction(t, k, prop; i=1)
    t, factor = tau_bose(t, prop.mesh[2].β)
    # return factor * Interp.interpND(view(prop.data, i, :, :), prop.mesh[2:3], (t, k))
    return factor * Interp.linear2D(view(prop.data, i, :, :), prop.mesh[2], prop.mesh[3], t, k)
end

function G0(t, k, param)
    β = param.β
    t, factor = tau_fermi(t, β)
    ε = k^2 / 2 / param.me - param.μ
    # result = green2(ε, t, β)
    result = factor * exp(-t * ε) / (exp(-ε * β) + 1)
    # if isnan(result) || isinf(result)
    #     println("t=$t, k=$k, ε=$ε")
    # end
    return result
end

function GG0(wn, k, param)
    β = param.β
    ε = k^2 / 2 / param.me - param.μ
    # ω = (2n + 1)π / β
    # return 1 / (ω^2 + ε^2)
    return 1 / (wn^2 + ε^2)
end

function calcF!(F, R0, Rt, param; norm=1.0)
    # compute F in imtime by dlr fourier transform
    Rdlr = to_dlr(Rt)
    Rw = dlr_to_imfreq(Rdlr)
    Fw = similar(Rw)
    for ind in eachindex(Fw)
        wn, k = Fw.mesh[1][ind[1]], Fw.mesh[2][ind[2]]
        Fw[ind] = real(Rw[ind] + R0[ind[2]]) * GG0(wn, k, param)
    end
    tgrid = F.mesh[1]
    Fdlr = to_dlr(Fw)
    Ft = dlr_to_imtime(Fdlr, tgrid)
    F.data .= Ft.data .* norm
    return true
end

function initR(param;
    mint=0.001,
    minK=0.001 * sqrt(param.T * param.me), maxK=10.0 * param.kF,
    order=6)

    # Euv = 100param.EF
    Euv = βEUV / param.β
    # rtol = 1e-10

    Nlog = floor(Int, nlog_factor * log10(param.β / mint))
    btgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, param.β], [0.0, param.β], Nlog, mint, order)
    tgrid = GreenFunc.ImTime(param.β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha, grid=btgrid)
    wn_mesh = GreenFunc.ImFreq(param.β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)

    Nk = floor(Int, 2.0 * log10(maxK / minK))
    kgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)

    ri = GreenFunc.MeshArray(kgrid; dtype=ComplexF64)

    R_freq = GreenFunc.MeshArray(wn_mesh, kgrid; dtype=Float64)
    for ind in eachindex(R_freq)
        Ω_c = 0.01
        ni, ki = ind[1], ind[2]
        ωn, k = wn_mesh[ni], kgrid[ki]
        e = k^2 / 2 / param.me - param.μ
        R_freq[ni, ki] = (1.0 - 2.0 * ωn^2 / (ωn^2 + Ω_c^2)) / (Ω_c^2 + e^2)
    end
    rdlr = to_dlr(R_freq)
    rtd = dlr_to_imtime(rdlr, tgrid)
    rt = dlr_to_imtime(rdlr)
    # rt = GreenFunc.MeshArray(tgrid, kgrid; dtype=ComplexF64)
    # ri[:] .= 0.0
    ri[:] .= 1.0
    # rt[:] .= 0.0

    return ri, rt, rtd
end

function loadR(fname, param;
    mint=0.001,
    minK=0.001 * sqrt(param.T * param.me), maxK=10.0 * param.kF,
    order=6)

    # Euv = 100param.EF
    Euv = βEUV / param.β
    # rtol = 1e-10

    Nlog = floor(Int, nlog_factor * log10(param.β / mint))
    btgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, param.β], [0.0, param.β], Nlog, mint, order)
    tdgrid = GreenFunc.ImTime(param.β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha, grid=btgrid)
    tgrid = GreenFunc.ImTime(param.β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)
    # wn_mesh = GreenFunc.ImFreq(param.β, FERMION; Euv=Euv, rtol=rtol, symmetry=:pha)

    println(tgrid)
    println(tdgrid)

    Nk = floor(Int, 2.0 * log10(maxK / minK))
    kgrid = Propagators.CompositeG.LogDensedGrid(:cheb, [0.0, maxK], [0.0, param.kF], Nk, minK, order)

    ri = GreenFunc.MeshArray(kgrid; dtype=ComplexF64)

    rt = GreenFunc.MeshArray(tgrid, kgrid; dtype=ComplexF64)
    rtd = GreenFunc.MeshArray(tdgrid, kgrid; dtype=ComplexF64)

    # load from file
    data = jldopen(fname, "r")
    R_ins = data["R_ins"]
    R_freq = data["R_freq"]

    rdlr = to_dlr(R_freq)
    Rt = dlr_to_imtime(rdlr, tgrid)
    Rtd = dlr_to_imtime(rdlr, tdgrid)

    for (ik, k) in enumerate(kgrid)
        ri[ik] = Interp.interp1D(view(R_ins.data, 1, :), R_ins.mesh[2], k)
        for (it, t) in enumerate(tgrid)
            rt[it, ik] = Interp.interp1D(view(Rt, it, :), Rt.mesh[2], k)
        end
        for (it, t) in enumerate(tdgrid)
            rtd[it, ik] = Interp.interp1D(view(Rtd, it, :), Rtd.mesh[2], k)
        end
    end
    println(rt[1, :])
    println(rtd[1, :])
    println("R0=$(real(R0(ri, rt, param)))")
    return ri, rt, rtd
end

function save_pcf(fname, param, ri, rt; R0=nothing, result=nothing)
    jldopen(fname, "w") do file
        file["param"] = param
        file["R_ins"] = ri
        file["R_imtime"] = rt
        if !isnothing(R0)
            file["R0"] = R0
        end
        if !isnothing(result)
            file["result"] = result
        end
    end
    return true
end

function R0(ri, rt, param)
    rdlr = Propagators.to_dlr(rt)
    rw = Propagators.to_imfreq(rdlr)
    kF = param.kF
    kgrid = rw.mesh[2]
    ikF = searchsortedfirst(kgrid, kF)
    # return 1.0 .+ ri[ikF] .+ rw[:, ikF]
    return ri[ikF] .+ rw[:, ikF]
end

function response(k, ri; norm=1)
    # return 1.0 + Interp.interp1D(view(ri.data, :), ri.mesh[1], k) / norm
    # return 1.0 + Interp.linear1D(view(ri.data, :), ri.mesh[1], k) / norm
    return Interp.linear1D(view(ri.data, :), ri.mesh[1], k) / norm
end

function response(t, k, rt; norm=1)
    t, factor = tau_fermi(t, rt.mesh[1].β)
    # return factor * Interp.interpND(view(rt.data, :, :), rt.mesh[:], (t, k)) / norm
    return factor * Interp.linear2D(view(rt.data, :, :), rt.mesh[1], rt.mesh[2], t, k) / norm
end

function diagram_gen(rs, beta)
    paramc = UEG.ParaMC(rs=rs, beta=beta)
    kF = paramc.kF
    partition = [(2, 0, 0)]
    channel = [PHr,]
    filter = [NoHartree, NoBubble]
    diagram = Ver4.diagram(paramc, partition; channel=channel, filter=filter)
    return paramc, diagram
end

function diagram_eval(diagram, K, T, paramc; idx=1, Fs=1.0, Fd=1.0)
    partition, diagpara, diag, root, extT = diagram
    d = diag[idx]
    weight = d.node.current
    ExprTree.evalKT!(d, K, T, paramc)
    w = sum(weight[r] for (ri, r) in enumerate(d.root))
    return w
end

end

using Test, BenchmarkTools

if abspath(PROGRAM_FILE) == @__FILE__
    using .Propagators

    # @testset "Propagators" begin

    #     param = Propagators.Parameter.defaultUnit(0.01, 1.0)

    #     rpai, rpat = Propagators.rpa(param)

    #     println(size(rpat))
    #     println(size(rpai))

    #     p = (0.00372, 0.0733)
    #     println(Propagators.interaction(p..., rpat))
    #     println(Propagators.interaction(p[2], rpai))
    #     @time Propagators.interaction(p..., rpat)
    #     @time Propagators.interaction(p[2], rpai)
    #     @time Propagators.G0(p..., param)

    #     @time Propagators.interaction(p..., rpat)
    #     @time Propagators.interaction(p[2], rpai)
    #     @time Propagators.G0(p..., param)

    #     Ri, Rt = Propagators.initR(param)
    #     @time Propagators.response(p..., Rt)
    #     @time Propagators.response(p[2], Ri)

    #     @time Propagators.response(p..., Rt)
    #     @time Propagators.response(p[2], Ri)

    #     pht = Propagators.phonon(param)
    #     @time Propagators.interaction(p..., pht)
    #     @time Propagators.interaction(p..., pht)
    # end

    # @testset "Test ElectronLiquid" begin
    #     using .Propagators.ElectronLiquid
    #     θ, rs = 0.1, 3.0
    #     mass2 = 1e-9
    #     Fs = -0.0
    #     beta = 1 / θ
    #     paramc = UEG.ParaMC(rs=rs, beta=beta,
    #         Fs=Fs, order=2,
    #         mass2=mass2, isDynamic=true)
    #     param = paramc.basic
    #     rpai, rpat = Propagators.rpa(param)
    #     p = (5.0, 1.0)
    #     V = 1 / real(Propagators.interaction(p[2], rpai))
    #     W = real(Propagators.interaction(p..., rpat) * V)
    #     println("V=$V, W=$W")

    #     V0 = UEG.Coulombinstant(p[2], paramc)
    #     println(V0)
    #     println(UEG.interactionDynamic(paramc, p[2], 0.0, p[1]))
    #     println(UEG.interactionStatic(paramc, p[2], 0.0, p[1]))
    # end

    # @testset "init and load R" begin
    #     fname = "run/data/PCFdata_5008.jld2"
    #     θ, rs = 0.002, 0.5
    #     param = Propagators.Parameter.rydbergUnit(θ, rs, 3)
    #     mint = 0.0001 * param.β
    #     minK, maxK = 0.1 * sqrt(param.T * param.me), 10param.kF
    #     order = 3
    #     Ri, Rt, Rtd = Propagators.initR(param; mint=mint, minK=minK, maxK=maxK, order=order)
    #     println(size(Ri))
    #     println(size(Rt))
    #     println(size(Rtd))
    #     Ri, Rt, Rtd = Propagators.loadR(fname, param; mint=mint, minK=minK, maxK=maxK, order=order)
    #     println(size(Ri))
    #     println(size(Rt))
    #     println(size(Rtd))
    #     println(Rtd.data[1, 1])

    #     calcF!(Rtd, Ri, Rt, param)
    #     println(Rtd.data[1, 1])

    # end

    @testset "diagram gen" begin
        θ, rs = 0.02, 0.5
        paramc, diagram = diagram_gen(rs, 1 / θ)
    end
end