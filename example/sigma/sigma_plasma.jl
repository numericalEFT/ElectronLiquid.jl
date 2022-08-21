using ElectronGas
using ElectronLiquid
using LinearAlgebra
using Measurements
using CompositeGrids
using Plots
using LaTeXStrings

const rs = 5.0
const beta = 25.0
const mass2 = 0.01

@inline function UEG.polarKW(q, n::Int, para::ParaMC)
    # return Polarization.Polarization0_ZeroTemp(q, n, para.basic) * para.spin * para.massratio
    # return Polarization.Polarization0_3dZeroTemp_LinearDispersion(q, n, para.basic) * para.spin * para.massratio
    if q < 1e-6
        q = 1e-6
    end
    wn = (2π * n) / para.β
    vF = para.kF / para.me
    x = abs(wn / (vF * q))
    factor = para.additional[1]
    return (1.0 / (factor * x^2 + 1.0)) * (-para.NF)
    # Π = density * (x^2.0 / (3 + x^2.0))
end

function Sw(val, kgrid, para)
    sw(d) = imag(d[2] - d[1]) / (2π / para.β)
    return [sw(val[:, ki]) for ki in eachindex(kgrid.grid)]
end

function Sk(val, kgrid, para)
    return [real(val[1, ki] - val[1, 1]) / (k^2 / (2 * para.me)) for (ki, k) in enumerate(kgrid.grid)]
end

const para = UEG.ParaMC(beta=25.0, rs=5.0, mass2=0.00001, isDynamic=true)
push!(para.additional, 0.0)
const kF = para.kF

const diag = Sigma.diagram(para, [(1, 0, 0),])

const factorList = [0.01,]

Nk, korder = 4, 2
minK = 0.2 * para.kF
const qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 2.0kF], [0.0, kF,], Nk, minK, korder)
# const qgrid = para.qgrid

plasmafreq(factor) = round(1 / sqrt(factor / 3), sigdigits=3)

sw = zeros(Measurement{Float64}, length(factorList), length(qgrid))
sk = zeros(Measurement{Float64}, length(factorList), length(qgrid))

for fi in eachindex(factorList)
    para0 = UEG.ParaMC(beta=beta, rs=rs, mass2=mass2, isDynamic=true)
    push!(para0.additional, factorList[fi])
    data, result = Sigma.KW(para0, diag; neval=1e7, kgrid=qgrid, ngrid=[-1, 0])
    if isnothing(data) == false
        val = data[(1, 0, 0)]
        sw[fi, :] .= Sw(val, qgrid, para0)
        sk[fi, :] .= Sk(val, qgrid, para0)
    end
end

let
    p = plot(xlabel=L"$k/k_F$", ylabel=L"$\frac{d\Sigma(k, w=0)}{d\omega}$", legend=:topright)
    for (fi, factor) in enumerate(factorList)
        plot!(p, qgrid ./ para.kF, sw[fi, :], label=L"$\tilde{\omega_p} = %$(plasmafreq(factor))\omega_p$")
    end
    savefig(p, "sigma_locality_z.pdf")
    display(p)
end

let
    p = plot(xlabel=L"$k/k_F$", ylabel=L"$\operatorname{Re}\frac{\Sigma(k, \omega_0)- \Sigma(0, \omega_0)}{k^2/2m}$", legend=:bottomright, xlim=[0.2, qgrid[end] / para.kF], ylim=[-1.0, 0.0])
    for (fi, factor) in enumerate(factorList)
        plot!(p, qgrid ./ para.kF, sk[fi, :], label=L"$\tilde{\omega_p} = %$(plasmafreq(factor))\omega_p$")
    end
    savefig(p, "sigma_locality_massratio.pdf")
    display(p)
end
