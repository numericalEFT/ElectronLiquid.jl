# include("../common/counterterm.jl")
# include("../common/para_builder.jl")
# using .UEG
# using .CounterTerm
using ElectronLiquid

para = UEG.ParaMC(
    rs=1.0,
    # Fs=-0.585,
    Fs=-0.0,
    beta=25.0,
    mass2=0.001,
    order=3,
    isDynamic=true,
)
Zrenorm = false
# Zrenorm = true

df = CounterTerm.fromFile()
mu, sw = CounterTerm.getSigma(df, UEG.paraid(para), para.order)
dmu, dz = CounterTerm.derive_onebody_parameter_from_sigma(para.order, mu, sw, zrenorm=Zrenorm)
z = CounterTerm.chemicalpotential_renormalization(para.order, sw, dmu)
# z[2] *= -1
# z[3] *= -1
println("dz = ", z)
sumz = accumulate(+, z)
if Zrenorm == false
    z = @. 1.0 / (1.0 + sumz)
else
    z = @. 1.0 - sumz
end
println("z = ", z)