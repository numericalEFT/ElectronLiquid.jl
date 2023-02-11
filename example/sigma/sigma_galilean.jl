using ElectronGas
using GreenFunc
using Plots
using LaTeXStrings

para = Parameter.rydbergUnit(1.0 / 25.0, 5.0)

sigma = SelfEnergy.G0W0(para)

sigma_w = sigma[1] |> to_dlr |> to_imfreq

kgrid = sigma_w.mesh[2]
ngrid = sigma_w.mesh[1]

z = SelfEnergy.zfactor(para, sigma[1])[1]
factor = 1 / z - 1

# kFidx = searchsortedfirst(kgrid, para.kF)
# plot(ngrid, real.(sigma_w[:, 1]), xlims=(-para.EF * 5, para.EF * 5), label="G0W0")
# plot(ngrid, imag.(sigma_w[:, kFidx]), xlims=(-para.EF * 5, para.EF * 5), label="G0W0")
# p = plot(xlims=(-para.EF * 5, para.EF * 5), ylims=[-0.075, 0.075])

sigma_w2 = deepcopy(sigma_w)
for kidx in eachindex(kgrid)
    # kidx = searchsortedfirst(kgrid, k)
    # isw = imag.(sigma_w[:, kidx])
    # rsw = real.(sigma_w[:, kidx] .+ sigma[1][kidx]) .+ factor * kgrid[kidx]^2 / 2 / para.me
    sigma_w2 .+= real.(sigma[1][kidx])
    # plot!(p, ngrid / para.EF, imag.(sigma_w[:, kidx]), label=L"$k/k_F=%$(k/para.kF)$")
end

########### imaginary part ################################
p = plot(xlims=(-para.EF * 5, para.EF * 5), ylims=[-0.075, 0.075], xlabel=L"$i\omega_n$", ylabel=L"$\mathrm{Im}\Sigma(k, i\omega_n)$")
for k in [0.0, para.kF / 3, para.kF / 3 * 2, para.kF, para.kF * 1.5, para.kF * 2]
    kidx = searchsortedfirst(kgrid, k)
    plot!(p, ngrid / para.EF, imag.(sigma_w[:, kidx]), label=L"$k/k_F=%$(k/para.kF)$")
end
display(p)

########## real part #################################
# p = plot(xlims=(-para.EF * 5, para.EF * 5), xlabel=L"$i\omega_n$", ylabel=L"$\mathrm{Re}\Sigma(k, i\omega_n)$")
# for k in [0.0, para.kF / 3, para.kF / 3 * 2, para.kF, para.kF * 1.5, para.kF * 2]
#     kidx = searchsortedfirst(kgrid, k)
#     plot!(p, ngrid / para.EF, real.(sigma_w[:, kidx] .+ sigma[1][kidx]) .+ factor * kgrid[kidx]^2 / 2 / para.me, label="k/kF=$(k/para.kF)")
# end
# display(p)

# ####################### difference ####################################
# new_kgrid = kgrid.grid .^ 2 / 2 / para.me

# p = plot(xlims=(0, para.EF * 2), ylims=[0.0, para.EF * 2], aspect_ratio=1.0, ylabel=L"$i\omega_n$", xlabel=L"$k^2/2m$")
# Plots.contour!(p, new_kgrid, ngrid, imag.(sigma_w2))
# display(p)