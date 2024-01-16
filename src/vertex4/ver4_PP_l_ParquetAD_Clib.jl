
function MC_PP_ParquetAD_Clib(
    para; kamp=[para.kF,], kamp2=kamp, n=[[0, 1, -1],], l=0,
    neval=1e6, filename::Union{String,Nothing}=nothing, reweight_goal=nothing,
    filter=[NoHartree, NoBubble],
    channel=[PHr, PHEr, PPr],
    partition=UEG.partition(para.order),
    print=0
)
    dim, β, kF, order = paras[1].para.dim, paras[1].para.β, paras[1].para.kF, paras[1].para.order
    Nw = length(paras[1].ωn)
    partition, diagpara, extT_labels, spin_conventions = diagram
    for p in paras
        # p.para.isDynamic && UEG.MCinitialize!(p.para)
        if p.para.isDynamic
            if NoBubble in diagpara[1].filter
                UEG.MCinitialize!(p.para, false)
            else
                UEG.MCinitialize!(p.para, true)
            end
            root_dir = joinpath(@__DIR__, "source_codeParquetAD/dynamic/")
        end
        @assert length(p.ωn) == Nw "All parameters must have the same frequency list"
        @assert p.para.dim == dim "All parameters must have the same dimension"
        @assert p.para.β ≈ β "All parameters must have the same inverse temperature"
        @assert p.para.kF ≈ kF "All parameters must have the same Fermi momentum"
    end

    @assert length(diagpara) == length(extT_labels) == length(spin_conventions)

    MaxLoopNum = maximum([key[1] for key in partition]) + 2

    df = CSV.read(root_dir * "loopBasis_ParquetADmaxOrder$(order).csv", DataFrame)
    loopBasis = [df[!, col] for col in names(df)]
    momLoopPool = FrontEnds.LoopPool(:K, dim, loopBasis)

    if paras[1].para.isDynamic
        f = jldopen(root_dir * "leafinfo_O$(order).jld2", "r")
        leafstates = f["leafstates"]
        leafvalues = f["values"]
        for l in leafvalues
            l = l .* 0.0
        end
    else
        leafstates = Vector{Vector{Ver4.LeafStateADVer4Dynamic}}()
        leafvalues = Vector{Vector{Float64}}()
        for key in partition
            key_str = join(string.(key))
            df = CSV.read(root_dir * "leafinfo_O$(order)_Parquet$key_str.csv", DataFrame)
            leafstates_par = Vector{Ver4.LeafStateADVer4Dynamic}()
            for row in eachrow(df)
                push!(leafstates_par, Ver4.LeafStateADVer4Dynamic(row[2], Sigma._StringtoIntVector(row[3]), row[4:end]..., 1))
            end
            push!(leafstates, leafstates_par)
            push!(leafvalues, df[!, names(df)[1]])
        end
    end

    root = zeros(Float64, maximum(length.(extT_labels)))
    # println(root)

    # all momentum will be sampled around the dimensionless Fermi momentum 1.0
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=3) # the first three momenta are external
    K.data[:, 1] = UEG.getK(kF, dim, 1)
    K.data[:, 2] = UEG.getK(kF, dim, 1)
    K.data[:, 3] = UEG.getK(kF, dim, 1)
    # all time variables will be sampled within [0.0, 1.0]
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha) # the first one is external
    T.data[1] = 0.0
    if dim == 3
        X = MCIntegration.Continuous(-1.0, 1.0, alpha=alpha) #x=cos(θ)
    elseif dim == 2
        X = MCIntegration.Continuous(0.0, 2π, alpha=alpha) #x=θ
    end
    N = MCIntegration.Discrete(1, length(paras), alpha=alpha) #index of paras

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = [zeros(ComplexF64, 2, Nw, length(paras)) for p in diagpara]

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = MCIntegration.Configuration(;
            var=(K, T, X, N),
            dof=dof,
            obs=obs,
            type=Weight,
            measurefreq=measurefreq,
            userdata=(paras, MaxLoopNum, extT_labels,
                spin_conventions, leafstates, leafvalues, momLoopPool,
                root, isLayered2D, partition),
            kwargs...
        )
    end
    result = integrate(integrand_ParquetAD_Clib; measure=measure_ParquetAD_Clib, config=config, solver=:mcmc, neval=neval, print=print, kwargs...)

    # function info(idx, di)
    #     return @sprintf("   %8.4f ±%8.4f", avg[idx, di], std[idx, di])
    # end

    if isnothing(result) == false
        if print >= 0
            report(result.config)
            # report(result; pick=o -> (real(o[1, 1, 1])), name="uu")
            # report(result; pick=o -> (real(o[2, 1, 1])), name="ud")
        end

        datadict = Dict{eltype(partition),Any}()
        for k in 1:length(dof)
            avg, std = result.mean[k], result.stdev[k]
            r = measurement.(real(avg), real(std))
            i = measurement.(imag(avg), imag(std))
            data = Complex.(r, i)
            datadict[partition[k]] = data
        end
        return datadict, result
    else
        return nothing, nothing
    end
end
