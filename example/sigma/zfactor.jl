# include("../common/counterterm.jl")
# include("../common/para_builder.jl")
# using .UEG
# using .CounterTerm
using ElectronLiquid
using Printf

rs = [1.0, 2.0, 3.0, 5.0, 6.0, 8.0]
mass2 = [0.01, 0.001]
Fs = [-0.0,]
beta = [25.0, ]
order = [3,]

const filename = "para.csv"
# const filename = "para_wn_1minus0.csv"

cache =[]

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=true)

    # Zrenorm = false
    # Zrenorm = true

    mu, sw = CounterTerm.getSigma(para, parafile=filename)
    dzi, dmu, dz = CounterTerm.sigmaCT(para.order, mu, sw)

    println("δz = ", dz)
    println("δz_inverse = ", dzi)
    sumzi = accumulate(+, dzi)
    # if Zrenorm == false

    z1 = @. 1.0 / (1.0 + sumzi)
    # z1 = sumzi
    # println("z without renormalization = $z")
    push!(cache, sumzi[3])

    # else
    sumz = accumulate(+, dz)
    z2 = @. 1.0 + sumz
    # end
    # println("z with renormalization = $z")
    printstyled("$(UEG.short(para))\n", color=:green)
    printstyled(@sprintf("%8s   %24s    %24s     %24s\n", "order", "chemical-potential", "no-z-renorm", "z-renorm"), color=:yellow)
    for o in 1:para.order
        @printf("%8d   %24s   %24s     %24s\n", o, "$(dmu[o])", "$(z1[o])", "$(z2[o])")
    end
end

b1=(cache[1]-cache[2])/(sqrt(0.01)-sqrt(0.001))
b2=(cache[2]-cache[3])/(sqrt(0.001)-sqrt(0.0001))
println(b1,", ", b2)
c1, c2 = cache[2]-b1*sqrt(0.001), cache[3]-b2*sqrt(0.0001)
println(c1,", ", c2)
println(1.0/(1+c1),", ", 1/(1+c2))