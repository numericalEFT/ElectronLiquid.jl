function diagPara(type, isDynamic::Bool, spin, order, filter, transferLoop=nothing)
    inter = [Interaction(ChargeCharge, isDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    if type == VacuumDiag
        innerLoopNum = order + 1
    else
        innerLoopNum = order
    end

    if Proper in filter
        @assert !isnothing(transferLoop) "transferLoop must be provided if Proper is in filter"
    end

    if isnothing(transferLoop)
        return Parquet.DiagPara(
            type=type,
            innerLoopNum=innerLoopNum,
            hasTau=true,
            spin=spin,
            interaction=inter,
            filter=filter,
        )
    else
        return Parquet.DiagPara(
            type=type,
            innerLoopNum=innerLoopNum,
            hasTau=true,
            spin=spin,
            interaction=inter,
            filter=filter,
            transferLoop=transferLoop
        )
    end
end

function diagram_parquet_noresponse(diagtype::Union{DiagramType,Symbol}, paramc::ParaMC, _partition::Vector{T},
    extra_deriv_orders::Vector{Int}=[], extra_dep_funcs::Vector{Function}=[];
    channels=[PHr, PHEr, PPr, Alli], filter=[NoHartree],
    transferLoop=nothing, extK=nothing, optimize_level=0
) where {T}
    if diagtype isa Symbol
        diagtype = _diagtype(diagtype)
    end
    diagpara = Vector{DiagPara}()

    FeynGraphs = diagdict_parquet(diagtype, _partition, extra_deriv_orders, extra_dep_funcs;
        filter=filter, channels=channels, spin=paramc.spin, isDynamic=paramc.isDynamic,
        transferLoop=transferLoop, extK=extK, optimize_level=optimize_level)
    extT_labels = Vector{Vector{Int}}[]
    partition = sort(collect(keys(FeynGraphs)))

    for p in partition
        push!(diagpara, diagPara(type, paramc.isDynamic, paramc.spin, p[1], filter, transferLoop))
        push!(extT_labels, [collect(g.properties.extT) for g in FeynGraphs[p]])
    end
    return (partition, diagpara, FeynGraphs, extT_labels)
end

function diagram_parquet_response(diagtype::Union{DiagramType,Symbol}, paramc::ParaMC, _partition::Vector{T},
    extra_deriv_orders::Vector{Int}=[], extra_dep_funcs::Vector{Function}=[];
    channels=[PHr, PHEr, PPr, Alli], filter=[NoHartree],
    transferLoop=nothing, extK=nothing, optimize_level=0
) where {T}
    if diagtype isa Symbol
        diagtype = _diagtype(diagtype)
    end
    diagpara = Vector{DiagPara}()

    FeynGraphs = diagdict_parquet(diagtype, _partition, extra_deriv_orders, extra_dep_funcs;
        filter=filter, channels=channels, spin=paramc.spin, isDynamic=paramc.isDynamic,
        transferLoop=transferLoop, extK=extK, optimize_level=optimize_level)
    extT_labels = Vector{Vector{Int}}[]
    spin_conventions = Vector{Response}[]
    partition = sort(collect(keys(FeynGraphs)))

    for p in partition
        push!(diagpara, diagPara(type, paramc.isDynamic, paramc.spin, p[1], filter, transferLoop))
        push!(extT_labels, [collect(g.properties.extT) for g in FeynGraphs[p]])
        push!(spin_conventions, [g.properties.response for g in FeynGraphs[p]])
    end
    return (partition, diagpara, FeynGraphs, extT_labels, spin_conventions)
end

function diagdict_parquet(diagtype::Union{DiagramType,Symbol}, _partition::Vector{T},
    extra_deriv_orders::Vector{Int}=[], extra_dep_funcs::Vector{Function}=[];
    channels=[PHr, PHEr, PPr, Alli], filter=[NoHartree], spin=2.0,
    isDynamic=false, transferLoop=nothing, extK=nothing, optimize_level=0
) where {T}
    if diagtype isa Symbol
        diagtype = _diagtype(diagtype)
    end
    MaxOrder = maximum([p[1] for p in _partition])
    MinOrder = minimum([p[1] for p in _partition])
    dict_graphs = Dict{Tuple{Int,Int,Int},Vector{Graph}}()
    leaf_dep_funcs = [pr -> pr isa FrontEnds.BareGreenId, pr -> pr isa FrontEnds.BareInteractionId]
    append!(leaf_dep_funcs, extra_dep_funcs)

    for order in MinOrder:MaxOrder
        partition_ind = findall(p -> p[1] == order, _partition)
        partition_order = [_partition[i] for i in partition_ind]
        Max_GD_o = maximum([p[2] for p in partition_order])
        Max_ID_o = maximum([p[3] for p in partition_order])
        para = diagPara(diagtype, isDynamic, spin, order, filter, transferLoop)
        ver4df = Parquet.build(para, extK; channels=channels)
        optimize!(ver4df.diagram, level=optimize_level)
        renormalization_orders = [Max_GD_o, Max_ID_o, extra_deriv_orders...]
        dict_graph_order = taylorAD(ver4df.diagram, renormalization_orders, leaf_dep_funcs)
        for key in keys(ver4_dict)
            p = (order, key[1], key[2])
            if p in _partition
                dict_graphs[p] = dict_graph_order[key]
            end
        end
    end
    return dict_graphs
end

function _diagtype(type::Symbol)
    if type == :freeEnergy
        return Parquet.VacuumDiag
    elseif type == :sigma
        return Parquet.SigmaDiag
    elseif type == :green
        return Parquet.GreenDiag
    elseif type == :chargePolar
        return Parquet.PolarDiag
    elseif type == :vertex3
        return Parquet.Ver3Diag
    elseif type in [:vertex4, :vertex4I]
        return Parquet.Ver4Diag
    else
        error("$type is not implemented")
    end
end