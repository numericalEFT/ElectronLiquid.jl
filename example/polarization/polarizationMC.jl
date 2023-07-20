using ElectronLiquid
using CompositeGrids
using FeynmanDiagram
using JLD2
# using MPI

# MPI.Init()
dim = 3
rs = [1.0,]
mass2 = [1.0,]
Fs = [-0.0,]
beta = [40.0]
order = [4,]
neval = 1e10
isDynamic = false
response = ChargeCharge  # TODO: implement other responses

# Measuring instantaneous or time-dependent polarization?
measure_instantaneous = true
# measure_instantaneous = false

# number of uniform tgrid points to use when measure_instantaneous is false
n_tau = 1001

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
    kF = para.kF

    Nk, korder = 4, 4
    minK = 0.2kF
    kgrid =
        CompositeGrid.LogDensedGrid(
            :uniform,
            [0.0, 3.0kF],
            [0.0, 2.0kF],
            Nk,
            minK,
            korder,
        ).grid
    tgrid = measure_instantaneous ? [para.β - 1e-8] : collect(LinRange(1e-8, para.β - 1e-8, n_tau))

    # NOTE: Partitions of type (1, n_μ, n_λ) with n_λ ≥ 1 are zero for the polariation.
    #       Here we generate them anyway, as they will be automatically ignored in the
    #       diagram generation step.
    partition = UEG.partition(_order) # (offset = 1, as for self-energy partitions)

    neighbor = UEG.neighbor(partition)
    @time diagram = Polarization.diagram(para, partition; response=response)
    valid_partition = diagram[1]  # may differ in size from partition

    reweight_goal = [1.0, 1.0, 1.0, 1.0,
        2.0, 2.0, 2.0, 4.0, 4.0, 8.0, 2.0, 2.0, 2.0,
        4.0, 4.0, 8.0, 4.0, 4.0, 8.0, 8.0, 2.0]

    # Pad reweight_goal if needed
    pad_value = 2.0
    reweight_pad = repeat([pad_value], max(0, length(valid_partition) - length(reweight_goal) + 1))
    reweight_goal = [reweight_goal; reweight_pad]
    @assert length(reweight_goal) ≥ length(valid_partition) + 1

    polarization, result = Polarization.KT(para, diagram;
        neighbor=neighbor, reweight_goal=reweight_goal[1:(length(valid_partition)+1)],
        kgrid=kgrid, tgrid=tgrid, neval=neval, parallel=:thread)

    if isnothing(polarization) == false
        if measure_instantaneous
            savename = "data_instantaneous_polarization_$(response)_K.jld2"
        else
            savename = "data_polarization_$(response)_KT.jld2"
        end
        jldopen(savename, "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (tgrid, kgrid, polarization)
        end
    end
end

