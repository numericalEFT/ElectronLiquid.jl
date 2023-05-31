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

println(df_re)
pl = df_re[!, 1]
pl = [eval(Meta.parse(pl[i])) for i in eachindex(pl)]

_mu = Dict()
for (i, p) in enumerate(pl)
    _mu[p] = measurement(df_re[i, 2])
    if(p == (1,0,2) || p == (1,2,0)) 
        _mu[p] = _mu[p]/2
    end
end

# println(_mu)

_ms = Dict() 
for (i, p) in enumerate(pl)
    dm = (measurement(df_re[i, 5]) - measurement(df_re[i, 2])) / (para.kF*0.01) * para.me / para.kF
    _ms[p] = dm
    if (p == (1, 0, 2) || p == (1, 2, 0))
        _ms[p] = _ms[p] / 2
    end
end

for p in sort([k for k in pl])
    println("$p: μ = $(_mu[p]), m = $(_ms[p])")
end


dmi_t, dμi_t, dm_t = CounterTerm.sigmaCT(para.order, _mu, _ms)
println(dmi_t)

