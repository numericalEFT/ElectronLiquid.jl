include("propagators.jl")
using .Propagators
using .Propagators.ElectronLiquid
using .Propagators.ElectronLiquid.ElectronGas.CompositeGrids
using .Propagators.ElectronLiquid.ElectronGas.GreenFunc

using Plots
using LaTeXStrings
using JLD2
using Roots

err_factor = 1 / sqrt(5)

function percent_error(result)
    funcs = result.config.userdata
    param = funcs.param
    ri, rt = funcs.Ri, funcs.Rt
    e1, e2, e3, e4 = real(result.stdev)
    erri = sqrt.(e1 .^ 2 + e2 .^ 2)
    errt = sqrt.(e3 .^ 2 + e4 .^ 2)
    kF = param.kF
    kgrid = rt.mesh[2]
    ikF = searchsortedfirst(kgrid, kF)
    println(real(ri[ikF]), " ", erri[ikF])
    println(real(rt[:, ikF]))
    println(errt[:, ikF])
    Ei = erri[ikF] / abs(real(ri[ikF]))
    Et = sum(errt[:, ikF]) / sum(abs.(real(rt[:, ikF])))
    return Et
end

function load_mcpcf_data(uid, dir)
    fname = dir * "mcpcf_$uid.jld2"
    jldopen(fname, "r") do file
        param = file["param"]
        ri = file["R_ins"]
        rt = file["R_imtime"]
        R0 = file["R0"]
        result = file["result"]
        err = percent_error(result)
        println("rs=$(param.rs), beta=$(param.beta), err=$(err)")
        return param, ri, rt, R0, result, err
    end
end

function load_pcf_data(uid, dir)
    fname = dir * "PCFdata_$uid.jld2"
    jldopen(fname, "r") do file
        ri = file["R_ins"]
        rw = file["R_freq"]
        return ri, rw
    end
end

function load_mcpcf_list(uids, dir)
    params, ris, rts, R0s, results = [], [], [], [], []
    R00s = []
    errs = []
    for uid in uids
        param, ri, rt, R0, result, err = load_mcpcf_data(uid, dir)
        push!(params, param)
        push!(ris, ri)
        push!(rts, rt)
        push!(R0s, R0)
        push!(results, result)
        ri0, rw0 = load_pcf_data(uid, dir)
        kF = param.kF
        kgrid = rw0.mesh[2]
        ikF = searchsortedfirst(kgrid, kF)
        R00 = real(ri0[ikF] .+ rw0[:, ikF])
        push!(R00s, R00)
        push!(errs, err)
    end
    return params, ris, rts, R0s, R00s, results, errs
end

function extract_Rkf(param, R_ins, R_freq)
    ωgrid = R_freq.mesh[1]
    kgrid = R_freq.mesh[2]

    kF = param.kF
    ikF = searchsortedfirst(kgrid, kF)

    Rkf_ins = R_ins[ikF]
    Rkf_freq = R_freq[:, ikF]

    Rkf = real(Rkf_freq) .+ Rkf_ins
    return ωgrid, Rkf
end

function Nfind_zero(data, grid)
    wgrid = CompositeGrids.SimpleG.Arbitrary([w for w in grid])
    f(x) = CompositeGrids.Interp.interp1D(data, wgrid, x)
    w0 = find_zero(f, (wgrid[1], wgrid[end]))
    return w0
end

function measure_chi(F_freq::GreenFunc.MeshArray)
    F_dlr = F_freq |> to_dlr
    F_ins = real(dlr_to_imtime(F_dlr, [F_freq.mesh[1].representation.β,])) * (-1)
    kgrid = F_ins.mesh[2]
    integrand = view(F_ins, 1, :)
    return real(CompositeGrids.Interp.integrate1D(integrand, kgrid))
end

function diff_rss()
    dir = "run/data/"
    uids = [106, 206, 506, 1006, 3006, 5006]
    params, ris, rts, r0s, r00s, results, errs = load_mcpcf_list(uids, dir)
    errs = errs .* err_factor
    rss = [param.rs for param in params]
    R0s = [r0[1] for r0 in r0s]
    R00s = [r00[1] for r00 in r00s]
    println(rss)
    println(R0s)
    println(R00s)

    p = plot()
    ylims!(0.0, maximum(1 ./ R00s))
    plot!(p, rss, 1 ./ R0s, yerr=errs ./ R0s, label="RPA+cross")
    plot!(p, rss, 1 ./ R00s, label="RPA")
    xlabel!(p, L"$r_s$")
    ylabel!(p, L"$1/R_0$")
    display(p)
    readline()
    savefig(p, "run/diff_rss.pdf")
end

function diff_temps()
    dir = "run/data/"
    uids = [3003, 3004, 3005, 3006, 3007, 3008]
    params, ris, rts, r0s, r00s, results, errs = load_mcpcf_list(uids, dir)
    errs = errs .* err_factor
    temps = [param.T / param.EF for param in params]
    R0s = [r0[1] for r0 in r0s]
    R00s = [r00[1] for r00 in r00s]
    println(temps)
    println(R0s)
    println(R00s)

    p = plot(xscale=:log10)
    # ylims!(0.0, maximum(1 ./ R00s))
    plot!(p, temps, 1 ./ R0s, yerr=errs ./ R0s, label="RPA+cross")
    plot!(p, temps, 1 ./ R00s, label="RPA")
    xlabel!(p, L"$T/T_F$")
    ylabel!(p, L"$1/R_0$")
    display(p)
    readline()
    savefig(p, "run/diff_temps.pdf")
end

diff_rss()
diff_temps()