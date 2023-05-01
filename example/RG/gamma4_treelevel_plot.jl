using Plots
using LaTeXStrings
using JLD2
pgfplotsx()

include("gamma4_treelevel.jl")

rs = 8.0
mass2 = 1e-5
beta = 100.0

para = UEG.ParaMC(rs=rs, beta=beta, Fs=-0.0, order=1, mass2=mass2, isDynamic=true, isFock=false)
Λgrid = CompositeGrid.LogDensedGrid(:gauss, [1.0 * para.kF, 20 * para.kF], [para.kF,], 8, 0.01 * para.kF, 8)
println(length(Λgrid))
println(Λgrid)

fs, us, dfs, dus = gamma4_treelevel_RG(para, Λgrid; verbose=1)

jldopen("data_f.jld2", "a+") do f
    key = "$(UEG.short(para))"
    if haskey(f, key)
        @warn("replacing existing data for $key")
        delete!(f, key)
    end
    f[key] = (para, Λgrid, fs, us, dfs, dus)
end

KOgrid = collect(LinRange(1.0 * para.kF, 5.0 * para.kF, 100))

F_cs = Interp.interp1DGrid(fs, Λgrid, KOgrid)
U_cs = Interp.interp1DGrid(us, Λgrid, KOgrid)

F_ko, U_ko = gamma4_treelevel_KO(para, KOgrid; verbose=1)

p = plot(titlefontsize=18,
    guidefontsize=18,
    tickfontsize=18,
    legendfontsize=18,
    labelfontsize=18,
    ylabel="Quasiparticle Interaction", xlabel=L"\Lambda/k_F", legend=:topright, size=(800, 600),
    xlims=[0.0, 3.0],
    ylims=[0.0, 1.2],
    linewidth=2, thickness_scaling=1
)
# plot!(p, (Λgrid .- para.kF) ./ para.kF, linewidth=2, fs .* Λgrid, label=L"CS")
# plot!(p, (KOgrid .- para.kF) ./ para.kF, linewidth=2, F_ko .* KOgrid, label=L"KO")

############# plot f ##################
# plot!(p, (Λgrid .- para.kF) ./ para.kF, linewidth=2, fs, label=L"CS")
# plot!(p, (KOgrid .- para.kF) ./ para.kF, linewidth=2, F_ko, label=L"KO")

############ plot u ############
plot!(p, (KOgrid .- para.kF) ./ para.kF, linewidth=2, -F_cs, label=L"-F^{CS}_s")
plot!(p, (KOgrid .- para.kF) ./ para.kF, linewidth=2, U_cs, label=L"U^{CS}_s")
plot!(p, (KOgrid .- para.kF) ./ para.kF, linewidth=2, -F_ko, label=L"-F^{KO}_s=U^{KO}_s")
#plot a vertical line at 2kF
savefig(p, "FsUs.pdf")
display(p)
readline()