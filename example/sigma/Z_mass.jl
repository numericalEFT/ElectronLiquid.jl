using JLD2
using Measurements
# using Gaston
using PyPlot
using ElectronLiquid
using ElectronGas
using LaTeXStrings

const filenameZ = "data_2D/data_Z_beta50.0_lam1.0.jld2"
const filenameK = "data_2D/data_K_beta50.0_lam1.0.jld2"
# Vibrant qualitative colour scheme from https://personal.sron.nl/~pault/
cdict = Dict(["blue" => "#0077BB", "cyan" => "#33BBEE", "teal" => "#009988", "orange" => "#EE7733", "red" => "#CC3311", "magenta" => "#EE3377", "grey" => "#BBBBBB"]);

function FreqDerivative(data, β)
    # return @. (imag(data[3, 1]) - imagZ(data[2, 1])) / (2π / β)
    return @. (imag(data[2, 1]) - imag(data[1, 1])) / (2π / β)
end

function mu(data)
    return real(data[2, 1])
end

fz = jldopen(filenameZ, "r")
key = keys(fz)[1]
value = fz[key]
para, ngrid, kgrid, data = value[1], value[2], value[3], value[4]
printstyled(key, color=:yellow)
# printstyled(UEG.short(para), color=:yellow)
println()

for p in sort([k for k in keys(data)])
    # println("$p: μ = $(mu(data[p]))   z = $(FreqDerivative(data[p], para.β))")
    # println("$p: ImΣ0 = $(imag(data[p][2,1]))   ImΣ-1 = $(imag(data[p][1,1])) ImΣ1 = $(imag(data[p][3,1]))")
    println("$p: ImΣ0 = $(imag(data[p][2,1]))   ImΣ-1 = $(imag(data[p][1,1]))")
    # println("$p: ReΣ0 = $(real(data[p][2,1]))   ReΣ-1 = $(real(data[p][1,1])) ReΣ1 = $(real(data[p][3,1]))")

end

_mu = Dict()
for (p, val) in data
    _mu[p] = mu(val)
end
# println(_mu)
_z = Dict()
for (p, val) in data
    _z[p] = FreqDerivative(val, para.β)
end

dzi, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _z)
z_factor = zeros(typeof(dz[1]), para.order)

# println("dz: ", dzi)
δsn = dzi

for i in eachindex(dzi)
    z_factor[i] = 1 + sum(dz[1:i])
    # println("z[$i]: ", z_factor[i])
    println("δs_$i=", δsn[i])
end



fk = jldopen(filenameK, "r")
key = keys(fk)[1]
valuek = fk[key]
para, ngrid, kgrid, data = valuek[1], valuek[2], valuek[3], valuek[4]
dim, β, kF, me = para.dim, para.β, para.kF, para.me
kF_label = searchsortedfirst(kgrid, kF)
println(kF_label)


klength = Int((length(kgrid) - 1) / 2)
# klength = 3
println(klength)
_ms = [Dict() for i in 1:klength]
# _ms3 = Dict()

for (p, val) in data
    # push!(real.(dReSigma, val[1, :] .- val[1, 1]))
    ek2 = real.(val[1, :] .- val[1, kF_label]) ./ (kgrid .- kgrid[kF_label]) * para.me / para.kF
    # println("order: $p, mass2 = $((ek2[kF_label-1]+ek2[kF_label+1])/2)")
    println("order: $p")
    for δki in 0:klength
        if δki == 0
            println("δk= $(round((kgrid[kF_label+δki]-kF)/kF;digits=3))kF, Re(Σ_i(0,kF+δk))= $(real(val[1,kF_label+δki])),Im(Σ_i(0,kF+δk))= $(imag(val[1,kF_label+δki]))")
        else
            println("δk= -$(round((kgrid[kF_label+δki]-kF)/kF;digits=3))kF, Re(Σ_i(0,kF+δk))= $(real(val[1,kF_label-δki]))")
            println("δk= $(round((kgrid[kF_label+δki]-kF)/kF;digits=3))kF, Re(Σ_i(0,kF+δk))= $(real(val[1,kF_label+δki]))")
            _ms[δki][p] = (ek2[kF_label-δki] + ek2[kF_label+δki]) / 2
            # _ms[δki][p] = ek2[kF_label+δki]
        end
    end
end

δki = 1
println("δk= $(round((kgrid[kF_label+δki]-kF)/kF;digits=3))kF")
for p in sort([k for k in keys(data)])
    println("$p: μ = $(_mu[p])   m = $(_ms[δki][p])")
end


dmi, dμi, dm = [], [], []


for i in 1:klength
    dmi_t, dμi_t, dm_t = CounterTerm.sigmaCT(para.order, _mu, _ms[i])
    # println(dmi_t)
    δmn = dmi_t
    # δmn = CounterTerm.renormalization(para.order, _ms[i], dmu, dzi, nbody=1, zrenorm=true)
    push!(dmi, dmi_t)
    push!(dμi, dμi_t)
    println("δk=$(round((kgrid[kF_label+i]-kF)/kF;digits=3))kF")
    # δMn = [δmn[1], (δmn[1])^2 + δmn[2], (δmn[1])^3 + 2δmn[1] * δmn[2] + δmn[3]]
    δMn = [δmn[1], (δmn[1])^2 + δmn[2]]
    # δrn = [δMn[1] + δsn[1], δMn[2] + δMn[1] * δsn[1] + δsn[2], δMn[3] + δMn[2] * δsn[1] + δMn[1] * δsn[2] + δsn[3]]
    δrn = [δMn[1] + δsn[1], δMn[2] + δMn[1] * δsn[1] + δsn[2]]
    println(dmi_t)
    for j in eachindex(δrn)
        # push!(m_eff, 1.0 / (1.0 + sum(dmi[1:i])*z_factor[i])) 
        println("order $j: δm_$j = $(δmn[j])  (m*/m)_$j =", 1.0 + sum(δrn[1:j]))
    end
end

dReSigma = []

######## Plot 
# interp = pyimport("scipy.interpolate")
style = PyPlot.matplotlib."style"
style.use(["science", "std-colors"])
color = [cdict["blue"], cdict["red"]]
#cmap = get_cmap("Paired")
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 16
rcParams["font.family"] = "Times New Roman"
figure(figsize=(4, 4))
dmp = []
for o in 1:para.order
    y = sum([z.val for z in dReSigma[j]] for j in 1:o)
    #y = [z.val for z in dReSigma[o]]
    e = [z.err for z in dReSigma[o]]
    #y ./= para.EF
    #e ./=para.EF
    y ./= (kgrid) .^ 2 / (2 * para.me)
    e ./= (kgrid) .^ 2 / (2 * para.me)
    # println(zk[o])
    errorbar(kgrid / kF, y, yerr=e, color=color[o], capsize=4, fmt="o", markerfacecolor="none", label="Order $o")
    push!(dmp, measurement(y[kF_label], e[kF_label]))
    # yfit = signal.savgol_filter(y, 5, 3)
    x = kgrid / kF
    y[1:3] .= y[kF_label]
    e[1:3] .= e[kF_label]
    spl = interp.UnivariateSpline(x, y, w=1.0 ./ e)
    yfit = spl(x)
    plot(x, yfit, color=color[o], linestyle="--")
    #println(y)
    #println(yfit)
end
#xlim([kgrid[1] / kF, kgrid[end] / kF])
xlim([kgrid[1] / kF, 2.2])
#ylim([-1.0, 0.0])
xlabel(L"$k/k_F$")

# if Zrenorm
#     ylabel(L"$z \cdot (\Sigma_{k, i\omega_0}-\Sigma_{0, i\omega_0})/(k^2/2m)$")
# else
ylabel(L"$(\Sigma_{k, i\omega_0}-\Sigma_{0, i\omega_0})/(k^2/2m)$")
# end
legend(loc=2)
#savefig("sigmaK_rs$(para.rs)_Fs$(para.Fs)_$(para.dim)d.pdf")
println(dmp)
