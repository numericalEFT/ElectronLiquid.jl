include("../common/counterterm.jl")
include("../common/para_builder.jl")
using .UEG
using .CounterTerm

para = ParaMC(
    rs=5.0,
    Fs=-0.0,
    beta=100.0,
    mass2=1e-5,
    order=2
)
Zrenorm = false

df = fromFile()
mu, sw = getSigma(df, UEG.paraid(para), para.order)
dmu, dz = derive_onebody_parameter_from_sigma(para.order, mu, sw)
sw = mergeInteraction(sw)
z = chemicalpotential_renormalization(para.order, sw, dmu)
println(z)