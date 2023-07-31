function sigma(paras,)

    partition = UEG.partition(_order)

    diagram_df = Sigma.diagram(para, partition; dR=true)
    sigma_df, result = Sigma.Generic(plist, diagram_df;
        neighbor=neighbor,
        neval=neval)

    diagram = Sigma.diagram(para, partition; dR=false)
    sigma, result = Sigma.Generic(plist, diagram;
        neighbor=neighbor, reweight_goal=reweight_goal[1:length(partition)+1],
        neval=neval)

    if isnothing(sigma) == false
        jldopen("data_Z_$scheme.jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (ngrid, Λgrid, sigma, sigma_df)
        end
    end
end