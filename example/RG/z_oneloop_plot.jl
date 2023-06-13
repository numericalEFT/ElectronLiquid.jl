# using Plots
# using LaTeXStrings
using ElectronGas
using ElectronLiquid
using CompositeGrids
using JLD2
# pgfplotsx()

rs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
# rs = [4.0,]
mass2 = [1e-5,]
beta = [100.0,]
order = [1,]
neval = 1e7

function zfactor(data, β)
    return @. (imag(data[2, :]) - imag(data[1, :])) / (2π / β)
end

z_RG, z_RPA = [], []

for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)

    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=1, mass2=_mass2, isDynamic=true, isFock=false)

    f = jldopen("data_f.jld2", "r")
    key = "$(UEG.short(para))"
    para, Λgrid, fs, us, dfs, dus = f[key]

    f = jldopen("data_Z.jld2", "r")
    key = "$(UEG.short(para))"
    para, ngrid, Λgrid, sigma = f[key]

    println(Λgrid)
    println(fs)
    println(dfs)

    dz_df = zfactor(sigma[(1, 0, 0)], para.β)
    # println(dz_df[1], " - ", dz_df[end])

    dz = Interp.integrate1D(dz_df .* dfs / para.NF, Λgrid)
    # dz = Interp.integrate1D(dfs, Λgrid)

    sigdyn, sigint = SelfEnergy.G0W0(para.basic, int_type=:ko_const, Fs=-fs[1])
    z = SelfEnergy.zfactor(para.basic, sigdyn, ngrid=[-1, 0])[1]
    println("rs=$_rs with fs=$(fs[1]): z=$z")
    ds_dw = 1 / z - 1
    println(dz, " vs ", ds_dw)
    _z = exp(-ds_dw - dz)
    # _z = exp(-ds_dw)

    sigdyn, sigint = SelfEnergy.G0W0(para.basic, int_type=:rpa)
    z = SelfEnergy.zfactor(para.basic, sigdyn, ngrid=[-1, 0])[1]
    push!(z_RPA, z)
    push!(z_RG, _z)
end
println(" rs    z_RG      z_RPA")
for (ri, _rs) in enumerate(rs)
    println("$(_rs)  $(z_RPA[ri]))  $(z_RG[ri])")
end

