using Plots
using LaTeXStrings
pgfplotsx()

include("z_oneloop.jl")
include("gamma4_treelevel.jl")

rs = 8.0
mass2 = 1e-5
beta = 100.0

para = UEG.ParaMC(rs=rs, beta=beta, Fs=-0.0, order=1, mass2=mass2, isDynamic=true, isFock=false)
Λgrid = CompositeGrid.LogDensedGrid(:gauss, [1.0 * para.kF, 20 * para.kF], [para.kF,], 4, 0.01 * para.kF, 4)

println(Λgrid)

# fs, us = gamma4_treelevel_RG(para, Λgrid; verbose=1)
# println(fs)
fs = zero(Λgrid.grid)
z = zfactor_oneloop_RG(para, Λgrid, fs; verbose=1)
println(z)

