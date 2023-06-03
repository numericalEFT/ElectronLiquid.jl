using JLD2
using Measurements
# using Gaston
# using PyPlot
using ElectronLiquid
using ElectronGas
using CSV
using DataFrames

rs = 1.0
mass2 = 1.0
Fs = -0.0
beta = 40.0
order = 3
para = UEG.ParaMC(rs=rs, beta=beta, Fs=Fs, order=order, mass2=mass2, isDynamic=false)

df_re = DataFrame(CSV.File("re_sigma_iw0_rs=1.0_beta_ef=40.0_mass2=1.0_neval=5.0e10_dK=0.001.csv"))
df_im = DataFrame(CSV.File("im_sigma_iw0_rs=1.0_beta_ef=40.0_mass2=1.0_neval=5.0e10_dK=0.001.csv"))
dw_im = DataFrame(CSV.File("im_sigma_iwn_kF_rs=1.0_beta_ef=40.0_mass2=1.0_neval=5.0e10.csv"))

println(df_re)
println(dw_im)

pl = df_re[!, 1]
pl = [eval(Meta.parse(pl[i])) for i in eachindex(pl)]

_mu = Dict()
for (i, p) in enumerate(pl)
    _mu[p] = measurement(df_re[i, 2]) / (factorial(p[2])*factorial(p[3]))
end

_z = Dict()
for (i, p) in enumerate(pl)
    dΣdw = (measurement(dw_im[i, 3]) - measurement(dw_im[i, 2])) / (2π/para.β)
    _z[p] = dΣdw / (factorial(p[2]) * factorial(p[3]))
end

# println(_mu)

_ms = Dict() 
for (i, p) in enumerate(pl)
    dm = (measurement(df_re[i, 5]) - measurement(df_re[i, 2])) / (para.kF*0.01) * para.me / para.kF
    _ms[p] = dm / (factorial(p[2]) * factorial(p[3]))
end

for p in sort([k for k in pl])
    println("$p: μ = $(_mu[p]), ∂Σ/∂k = $(_ms[p]), ∂Σ/∂ω = $(_z[p])")
end

dzi, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _z)

z_factor = zeros(typeof(dz[1]), para.order)
δsn = dzi

for i in eachindex(dzi)
    z_factor[i] = 1 + sum(dz[1:i])
    # println("z[$i]: ", z_factor[i])
    println("δs_$i=", δsn[i])
end



