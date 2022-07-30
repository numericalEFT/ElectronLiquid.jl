include("../common/counterterm.jl")
include("../common/para_builder.jl")
using .UEG
using .CounterTerm

para = ParaMC(
    rs=5.0,
    Fs=-0.585,
    beta=25.0,
    mass2=0.01,
    order=3
)
Zrenorm = false

df = fromFile()
mu, sw = getSigma(df, UEG.paraid(para), para.order)
dmu, dz = derive_onebody_parameter_from_sigma(para.order, mu, sw, zrenorm=Zrenorm)
z = chemicalpotential_renormalization(para.order, sw, dmu)
println("dz = ", z)
sumz = accumulate(+, z)
if Zrenorm == false
    z = @. 1.0 / (1.0 + sumz)
else
    z = @. 1.0 - sumz
end
println("z = ", z)