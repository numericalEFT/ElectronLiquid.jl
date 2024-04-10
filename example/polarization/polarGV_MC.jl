using ElectronLiquid, FeynmanDiagram
using CompositeGrids
using JLD2, Printf

dim = 2
rs = [0.5,]
mass2 = [4.0,]
Fs = [-0.0,]
beta = [50.0]
order = [3,]
neval = 1e7
isDynamic = false
isFock = false
diagGenerate = :GV
# diagGenerate = :Parquet
response = ChargeCharge
# response = SpinSpin


for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock)
    kF = para.kF
    ######### calculate K dependence #####################
    Nk, korder = 4, 4
    minK = 0.2kF
    # kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 2.2kF], [kF,], Nk, minK, korder).grid
    # kgrid = kF .+ [-0.1, -0.05, -0.03, -0.01, -0.005, -0.001, 0, 0.001, 0.005, 0.01, 0.03, 0.05, 0.1] * kF
    kgrid = kF .+ [-0.1, -0.05, 0, 0.05, 0.1] * kF
    ngrid = [0,]

    # partition = [(1, 0, 0), (2, 0, 1), (2, 1, 0), (3, 0, 0)]
    _partition = UEG.partition(_order)
    reweight_goal = Float64[]
    partition = Vector{eltype(_partition)}()
    for (order, sOrder, vOrder) in _partition
        order == 1 && vOrder > 0 && continue
        push!(reweight_goal, 2.0^(2order + sOrder + vOrder - 2))
        push!(partition, (order, sOrder, vOrder))
    end
    push!(reweight_goal, 4.0)
    neighbor = UEG.neighbor(partition)

    if diagGenerate == :GV
        @time diagram = Polarization.diagramGV(para, partition; response=response)
        polar, result = Polarization.GV(para, diagram;
            neighbor=neighbor, reweight_goal=reweight_goal,
            kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:thread)
    elseif diagGenerate == :Parquet
        @time diagram = Polarization.diagram(para, partition; response=response)
        polar, result = Polarization.KW(para, diagram;
            neighbor=neighbor, reweight_goal=reweight_goal,
            kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:thread)
    else
        error("unknown diagrams' generated type")
    end

    if isnothing(polar) == false
        jldopen("data_$(response)Polar_K.jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (ngrid, kgrid, polar)
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            @printf("%10s  %10s   %10s   %10s   %10s \n", "q/kF", "real(avg)", "err", "imag(avg)", "err")
            r, i = real(polar[key]), imag(polar[key])
            for (in, n) in enumerate(ngrid)
                println("n = $n")
                for (iq, q) in enumerate(kgrid)
                    @printf("%10.6f  %10.6f ± %10.6f   %10.6f ± %10.6f\n", q[1] / kF, r[in, iq].val, r[in, iq].err, i[in, iq].val, i[in, iq].err)
                end
            end
        end
    end
end

