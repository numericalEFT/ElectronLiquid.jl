using JLD2
using Measurements
# using Gaston
# using PyPlot
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

# function ek(val, kgrid, kF_label)
#     return real(val[1, kF_label] - val[1, 1]) / (kgrid[kF_label]^2 / (2me))
# end

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
z_factor = zeros(typeof(dz[1]),para.order)
# println("dz: ", dzi)
for i in eachindex(dzi)
    z_factor[i] = 1.0 / (1.0 + sum(dzi[1:i]))
    println("z[$i]: ", 1.0 / (1.0 + sum(dzi[1:i])))  
end
println(dz)
δsn = dzi
# δsn = CounterTerm.renormalization(para.order, _z, dmu, dzi; nbody=1, zrenorm=true)
println(δsn)

fk = jldopen(filenameK, "r")
key = keys(fk)[1]
valuek = fk[key]
para, ngrid, kgrid, data = valuek[1], valuek[2], valuek[3], valuek[4]
dim, β, kF, me = para.dim, para.β, para.kF, para.me
kF_label = searchsortedfirst(kgrid, kF)
println(kF_label)


klength = Int((length(kgrid) - 1) / 2)
println(klength)
_ms = [Dict() for i in 1:klength]
# _ms3 = Dict()
for (p, val) in data
    ek2 = real.(val[1, :] .- val[1, kF_label]) ./ (kgrid .- kgrid[kF_label]) * para.me / para.kF
    # println("order: $p, mass2 = $((ek2[kF_label-1]+ek2[kF_label+1])/2)")
    println("order: $p")
    for δki in 0: klength
        if δki == 0
            println("δk= $(round((kgrid[kF_label+δki]-kF)/kF;digits=3))kF, Re(Σ_i(0,kF+δk))= $(real(val[1,kF_label+δki])),Im(Σ_i(0,kF+δk))= $(imag(val[1,kF_label+δki]))")
        else
            println("δk= -$(round((kgrid[kF_label+δki]-kF)/kF;digits=3))kF, Re(Σ_i(0,kF+δk))= $(real(val[1,kF_label-δki]))")
            println("δk= $(round((kgrid[kF_label+δki]-kF)/kF;digits=3))kF, Re(Σ_i(0,kF+δk))= $(real(val[1,kF_label+δki]))")
            _ms[δki][p] = (ek2[kF_label-δki] + ek2[kF_label+δki]) / 2
            _ms[δki][p] =  ek2[kF_label+δki]
        end
    end
end
δki = 3
println("δk= $(round((kgrid[kF_label+δki]-kF)/kF;digits=3))kF")
for p in sort([k for k in keys(data)])
    println("$p: μ = $(_mu[p])   m = $(_ms[δki][p])")
end


dmi, dμi, dm = [ ], [ ], [ ]


for i in 1:klength
    dmi_t, dμi_t, dm_t = CounterTerm.sigmaCT(para.order, _mu, _ms[i])
    # println(dmi_t)
    δmn = dmi_t
    # δmn = CounterTerm.renormalization(para.order, _ms[i], dmu, dzi, nbody=1, zrenorm=true)
    push!(dmi,dmi_t)
    push!(dμi, dμi_t)
    println("δk=$(round((kgrid[kF_label+i]-kF)/kF;digits=3))kF")
    δMn = [δmn[1], (δmn[1])^2 + δmn[2], (δmn[1])^3 + 2δmn[1]*δmn[2]+ δmn[3]]
    # δMn = [-δmn[1], (δmn[1])^2 - δmn[2], -(δmn[1])^3 + 2δmn[1]*δmn[2]- δmn[3]]
    δrn = [δMn[1] + δsn[1], δMn[2] + δMn[1] * δsn[1] + δsn[2], δMn[3] + δMn[2] * δsn[1] + δMn[1]*δsn[2] + δsn[3]]
    # δrn = [δMn[1] - δsn[1], δMn[2] - δMn[1] * δsn[1] - δsn[2], δMn[3] - δMn[2] * δsn[1] - δMn[1] * δsn[2] - δsn[3]]
    println(dmi_t)
    for j in eachindex(δrn)
    # push!(m_eff, 1.0 / (1.0 + sum(dmi[1:i])*z_factor[i])) 
        println("(m*/m) order $j: ", 1.0 + sum(δrn[1:j]))
    end
end

