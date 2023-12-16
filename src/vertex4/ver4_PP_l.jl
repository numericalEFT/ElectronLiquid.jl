
function MC_PP(para; kamp=[para.kF,], kamp2=kamp, n=[[0, 1, -1],], l=0,
    neval=1e6, filename::Union{String,Nothing}=nothing, reweight_goal=nothing,
    filter=[NoHartree, NoBubble],
    channel=[PHr, PHEr, PPr],
    partition=UEG.partition(para.order),
    print=0
)

    kF = para.kF
    _order = para.order

    # partition = UEG.partition(_order)


    diagram = Ver4.diagram(para, partition; channel=channel, filter=filter)

    partition = diagram[1] # diagram like (1, 1, 0) is absent, so the partition will be modified
    neighbor = UEG.neighbor(partition)

    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            # push!(reweight_goal, 8.0^(order + vOrder - 1))
            push!(reweight_goal, 8.0^(order - 1))
        end
        push!(reweight_goal, 1.0)
        println(length(reweight_goal))
    end

    paras = [Ver4.OneAngleAveraged(para, [kamp[1], kamp2[1]], n, :PP, l),]
    ver4, result = Ver4.one_angle_averaged(paras, diagram;
        neval=neval, print=print,
        neighbor=neighbor,
        reweight_goal=reweight_goal
    )

    if isnothing(ver4) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (kamp, n, l, ver4)
            end
        end
    end
    return ver4, result
end

function MC_PP_ParquetAD(para; kamp=[para.kF,], kamp2=kamp, n=[[0, 1, -1],], l=0,
    neval=1e6, filename::Union{String,Nothing}=nothing, reweight_goal=nothing,
    filter=[NoHartree, NoBubble],
    channel=[PHr, PHEr, PPr],
    partition=UEG.partition(para.order),
    isClib=true,
    print=0
)

    kF = para.kF
    _order = para.order

    # partition = UEG.partition(_order)

    if isClib
        diagram = Ver4.diagramParquet_load(para, partition; filter=filter)
    else
        diagram = Ver4.diagramParquet(para, partition; channel=channel, filter=filter)
    end

    partition = diagram[1] # diagram like (1, 1, 0) is absent, so the partition will be modified
    neighbor = UEG.neighbor(partition)

    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            # push!(reweight_goal, 8.0^(order + vOrder - 1))
            push!(reweight_goal, 8.0^(order - 1))
        end
        push!(reweight_goal, 1.0)
        println(length(reweight_goal))
    end

    paras = [Ver4.OneAngleAveraged(para, [kamp[1], kamp2[1]], n, :PP, l),]
    if isClib
        ver4, result = Ver4.one_angle_averaged_ParquetAD_Clib(paras, diagram;
            neval=neval, print=print,
            neighbor=neighbor,
            reweight_goal=reweight_goal
        )
    else
        ver4, result = Ver4.one_angle_averaged_ParquetAD(paras, diagram;
            neval=neval, print=print,
            neighbor=neighbor,
            reweight_goal=reweight_goal
        )
    end
    if isnothing(ver4) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (kamp, n, l, ver4)
            end
        end
    end
    return ver4, result
end
