using JLD2
using Measurements
# using Gaston
using PyPlot
using ElectronLiquid
using ElectronGas

const filenameZ = "data_Z.jld2"
const filenameK = "data_K.jld2"

function zfactor(data, β)
    return @. (imag(data[2, 1]) - imag(data[1, 1])) / (2π / β)
end

function mu(data)
    return real(data[1, 1])
end

function ek(val, kgrid, kF_label)
    return real(val[1, kF_label] - val[1, 1]) / (kgrid[kF_label]^2 / (2me))
end

fz = jldopen(filenameZ, "r")
key = keys(fz)[1]
value = fz[key]
para, ngrid, kgrid, data = value[1], value[2], value[3], value[4]
printstyled(UEG.short(para), color=:yellow)
println()
for p in sort([k for k in keys(data)])
    println("$p: μ = $(mu(data[p]))   z = $(zfactor(data[p], para.β))")
end

_mu = Dict()
for (p, val) in data
    _mu[p] = mu(val)
end
# println(_mu)
_z = Dict()
for (p, val) in data
    _z[p] = zfactor(val, para.β)
end
dzi, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _z)
# z_factor = [ ]
# println("dz: ", dzi)
for i in eachindex(dzi)
    # push!(zfactor, 1.0 / (1.0 + sum(dzi[1:i])))
    println("z[$i]: ", 1.0 / (1.0 + sum(dzi[1:i])))
end

fk = jldopen(filenameK, "r")
key = keys(fk)[1]
valuek = fk[key]
para, ngrid, kgrid, data = valuek[1], valuek[2], valuek[3], valuek[4]
dim, β, kF, me = para.dim, para.β, para.kF, para.me
kF_label = searchsortedfirst(kgrid, kF)
println(kF_label)

# _ms = Dict()
# for (p, val) in data
#     ek = real.(val[1, :] .- val[1, 1]) ./ (kgrid .^ 2 / (2me))
#     println("order: $p, mass1 = $(ek[kF_label])")
#     # println(ek)
#     _ms[p] = ek[kF_label]
# end

# dmi, dμ, dm = CounterTerm.sigmaCT(para.order, _mu, _ms)

# m_eff = []

# for i in eachindex(dmi)
#     # push!(m_eff, 1.0 / (1.0 + sum(dmi[1:i])*z_factor[i])) 
#     println("m_eff[$i]: ", 1.0 /((1.0 - sum(dmi[1:i])) / (1.0 + sum(dzi[1:i]))))
# end

δk = 1
δk2 = 4

_ms2 = Dict()
_ms3 = Dict()
for (p, val) in data
    ek2 = real.(val[1, :] .- val[1, kF_label]) ./ (kgrid .- kgrid[kF_label]) * para.me / para.kF
    # println("order: $p, mass2 = $((ek2[kF_label-1]+ek2[kF_label+1])/2)")
    println("order: $p")
    for δki in -3:3
        println("δk= $(round((kgrid[kF_label+δki]-kF)/kF;digits=3))kF, Re(Σ_i(0,kF+δk))= $(real(val[1,kF_label+δki]))")
    end
    # println(ek2)
    _ms2[p] = (ek2[kF_label-δk] + ek2[kF_label+δk]) / 2
    _ms3[p] = (ek2[kF_label-δk2] + ek2[kF_label+δk2]) / 2
end

dmi2, dμ2, dm2 = CounterTerm.sigmaCT(para.order, _mu, _ms2)
dmi3, dμ3, dm3 = CounterTerm.sigmaCT(para.order, _mu, _ms3)


# m_eff = []
# println("δk-=$(round(-(kgrid[kF_label-δk]-kF)/kF;digits=3))kF")
println("δk=$(round((kgrid[kF_label+δk]-kF)/kF;digits=3))kF")
for i in eachindex(dmi2)
    # push!(m_eff, 1.0 / (1.0 + sum(dmi[1:i])*z_factor[i])) 
    println("m_eff[$i]: ", 1.0 / ((1.0 - sum(dmi2[1:i])) / (1.0 + sum(dzi[1:i]))))
end

println("δk=$(round((kgrid[kF_label+δk2]-kF)/kF;digits=3))kF")
for i in eachindex(dmi3)
    # push!(m_eff, 1.0 / (1.0 + sum(dmi[1:i])*z_factor[i])) 
    println("m_eff[$i]: ", 1.0 / ((1.0 - sum(dmi3[1:i])) / (1.0 + sum(dzi[1:i]))))
end

