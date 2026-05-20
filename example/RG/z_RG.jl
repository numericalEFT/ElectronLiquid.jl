using ElectronGas
using ElectronLiquid
using CompositeGrids
using GreenFunc
using JLD2

include("gamma4_treelevel.jl")

function zfactor(data, β, ngrid)
    if ngrid == [0, 1]
        return @. (imag(data[3, :]) - imag(data[2, :])) / (2π / β)
    elseif ngrid == [-1, 0]
        return @. (imag(data[2, :]) - imag(data[1, :])) / (2π / β)
    else
        error("ngrid = $ngrid not implemented")
    end
end

function z_flow(para, scheme=:KO; ngrid=[0, 1], verbose=1)
    filename = "data_Z_$(scheme).jld2"

    f = jldopen(filename, "r")
    key = "$(UEG.short(para))"
    _ngrid, Λgrid, sigma, sigma_df = f[key]

    @assert _ngrid == [-1, 0, 1]

    if scheme == :KO
        fs, us, dfs, dus = gamma4_treelevel_KO(para, Λgrid)
    else
        error("scheme not implemented")
    end
    # println(fs)
    # println(dfs)

    dz_df = zfactor(sigma_df[(1, 0, 0)], para.β, ngrid)
    dz = [Interp.integrate1D(dz_df .* dfs / para.NF, Λgrid, (Λgrid[li], Λgrid[end])) for li in eachindex(Λgrid)]

    # sigdyn, sigint = SelfEnergy.G0W0(para.basic, Λgrid; int_type=:ko_const, Fs=-fs[1])
    # Σ_freq = dlr_to_imfreq(to_dlr(sigdyn), [-1, 0, 1])
    # ds_dw = zfactor(Σ_freq, para.β, ngrid)
    ds_dw = zfactor(sigma[(1, 0, 0)], para.β, ngrid)

    z_RG = @. exp(-ds_dw - dz)
    # z_RG = @. exp(-ds_dw)
    println(ds_dw[1] + dz[1])
    println(ds_dw[1], ", ", dz[1], ", ", z_RG[1])

    z_RPA = @. 1 / (1 + ds_dw)

    if verbose > 0
        kF_label = searchsortedfirst(Λgrid, para.kF)
        println("rs = $(para.rs), kF = $(para.kF), kF_label = $kF_label with z_RPA = $(z_RPA[kF_label]) and z_RG = $(z_RG[kF_label])")
    end
    return Λgrid, z_RG, ds_dw, dz_df
end

if abspath(PROGRAM_FILE) == @__FILE__

    # para = UEG.ParaMC(rs=4.0, beta=25.0, Fs=-0.0, order=1, mass2=0.01, isDynamic=true, isFock=false)

    # rs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
    rs = [4.0,]
    mass2 = [0.001,]
    _Fs = [-0.0,]
    beta = [25.0,]
    order = [1,]

    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, _Fs, beta, order)
        para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=1, mass2=_mass2, isDynamic=true, isFock=false)
        z_flow(para, :KO; ngrid=[-1, 0])
    end
end