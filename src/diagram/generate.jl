function diagPara(type, isDynamic::Bool, order, spin=2, filter=[NoHartree], transferLoop=nothing)
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
    leaf_dep_funcs::Vector{Function}=[pr -> pr isa FrontEnds.BareGreenId, pr -> pr isa FrontEnds.BareInteractionId];
    filter=[NoHartree], extK=nothing, optimize_level=0
) where {T}
    if diagtype isa Symbol
        diagtype = _diagtype(diagtype)
    end
    diagpara = Vector{DiagPara}()

    FeynGraphs = diagdict_parquet(diagtype, _partition, leaf_dep_funcs;
        filter=filter, spin=paramc.spin, isDynamic=paramc.isDynamic,
        extK=extK, optimize_level=optimize_level)
    extT_labels = Vector{Vector{Int}}[]
    partition = sort(collect(keys(FeynGraphs)))

    for p in partition
        push!(diagpara, diagPara(diagtype, paramc.isDynamic, p[1], paramc.spin, filter))
        push!(extT_labels, [collect(g.properties.extT) for g in FeynGraphs[p]])
    end
    return (partition, diagpara, FeynGraphs, extT_labels)
end

function diagram_parquet_response(diagtype::Union{DiagramType,Symbol}, paramc::ParaMC, _partition::Vector{T},
    leaf_dep_funcs::Vector{Function}=[pr -> pr isa FrontEnds.BareGreenId, pr -> pr isa FrontEnds.BareInteractionId];
    channels=[PHr, PHEr, PPr, Alli], filter=[NoHartree],
    transferLoop=nothing, extK=nothing, optimize_level=0
) where {T}
    if diagtype isa Symbol
        diagtype = _diagtype(diagtype)
    end
    diagpara = Vector{DiagPara}()

    FeynGraphs = diagdict_parquet(diagtype, _partition, leaf_dep_funcs;
        filter=filter, channels=channels, spin=paramc.spin, isDynamic=paramc.isDynamic,
        transferLoop=transferLoop, extK=extK, optimize_level=optimize_level)
    extT_labels = Vector{Vector{Int}}[]
    spin_conventions = Vector{Response}[]
    partition = sort(collect(keys(FeynGraphs)))

    for p in partition
        push!(diagpara, diagPara(diagtype, paramc.isDynamic, p[1], paramc.spin, filter, transferLoop))
        push!(extT_labels, [collect(g.properties.extT) for g in FeynGraphs[p]])
        push!(spin_conventions, [g.properties.response for g in FeynGraphs[p]])
    end
    return (partition, diagpara, FeynGraphs, extT_labels, spin_conventions)
end

function diagdict_parquet(diagtype::Union{DiagramType,Symbol}, _partition::Vector{T},
    leaf_dep_funcs::Vector{Function}=[pr -> pr isa FrontEnds.BareGreenId, pr -> pr isa FrontEnds.BareInteractionId];
    # extra_deriv_orders::Vector{Int} = Int[], extra_dep_funcs::Vector{Function}=Function[];
    channels=[PHr, PHEr, PPr, Alli], filter=[NoHartree], spin=2.0, isDynamic=false,
    transferLoop=nothing, extK=nothing, optimize_level=0
) where {T}
    deriv_num = length(leaf_dep_funcs)
    @assert length(_partition[1]) == deriv_num + 1 "partition should have $deriv_num+1 entries"

    if diagtype isa Symbol
        diagtype = _diagtype(diagtype)
    end
    MaxOrder = maximum([p[1] for p in _partition])
    MinOrder = minimum([p[1] for p in _partition])
    dict_graphs = Dict{NTuple{deriv_num + 1,Int},Vector{Graph}}()
    # leaf_dep_funcs = [pr -> pr isa FrontEnds.BareGreenId, pr -> pr isa FrontEnds.BareInteractionId]
    # append!(leaf_dep_funcs, extra_dep_funcs)

    for order in MinOrder:MaxOrder
        partition_ind = findall(p -> p[1] == order, _partition)
        partition_order = [_partition[i] for i in partition_ind]

        # Max_GD_o = maximum([p[2] for p in partition_order])
        # Max_ID_o = maximum([p[3] for p in partition_order])
        para = diagPara(diagtype, isDynamic, order, spin, filter, transferLoop)
        graph_df = Parquet.build(para, extK; channels=channels)
        optimize!(graph_df.diagram, level=optimize_level)
        optimize!(graph_df.diagram, level=optimize_level)

        renormalization_orders = Int[]
        for i in 1:deriv_num
            push!(renormalization_orders, maximum([p[i+1] for p in partition_order]))
        end
        # renormalization_orders = [Max_GD_o, Max_ID_o, extra_deriv_orders...]
        dict_graph_order = taylorAD(graph_df.diagram, renormalization_orders, leaf_dep_funcs)
        for key in keys(dict_graph_order)
            p = (order, key...)
            if p in _partition
                dict_graphs[p] = dict_graph_order[key]
            end
        end
    end
    return dict_graphs
end

function diagram_GV_freeE(paramc::ParaMC, _partition::Vector{T},
    leaf_dep_funcs::Vector{Function}=[pr -> pr isa FrontEnds.BareGreenId, pr -> pr isa FrontEnds.BareInteractionId];
    filter=[NoHartree], optimize_level=0
) where {T}
    deriv_num = length(leaf_dep_funcs)
    @assert length(_partition[1]) == deriv_num + 1 "partition should have $deriv_num+1 entries"

    diagpara = Vector{DiagPara}()
    MaxOrder = maximum([p[1] for p in _partition])
    MinOrder = minimum([p[1] for p in _partition])
    dict_graphs = Dict{NTuple{deriv_num + 1,Int},Vector{Graph}}()
    spinPolarPara = 2.0 / paramc.spin - 1

    for order in MinOrder:MaxOrder
        partition_ind = findall(p -> p[1] == order, _partition)
        partition_order = [_partition[i] for i in partition_ind]

        diagrams = GV.diagsGV(:freeEnergy, order, spinPolarPara=spinPolarPara, filter=filter)
        optimize!(diagrams, level=optimize_level)
        optimize!(diagrams, level=optimize_level)

        renormalization_orders = Int[]
        for i in 1:deriv_num
            push!(renormalization_orders, maximum([p[i+1] for p in partition_order]))
        end
        # renormalization_orders = [Max_GD_o, Max_ID_o, extra_deriv_orders...]
        dict_graph_order = taylorAD(diagrams, renormalization_orders, leaf_dep_funcs)
        for key in keys(dict_graph_order)
            p = (order, key...)
            if p in _partition
                dict_graphs[p] = dict_graph_order[key]
            end
        end
    end

    partition = sort(collect(keys(dict_graphs)))
    for p in partition
        push!(diagpara, diagPara(VacuumDiag, paramc.isDynamic, p[1], paramc.spin, filter))
    end
    return (partition, diagpara, dict_graphs)
end

function diagram_GV_noresponse(diagtype::Symbol, paramc::ParaMC, _partition::Vector{T},
    leaf_dep_funcs::Vector{Function}=[pr -> pr isa FrontEnds.BareGreenId, pr -> pr isa FrontEnds.BareInteractionId];
    filter=[NoHartree], optimize_level=0
) where {T}
    deriv_num = length(leaf_dep_funcs)
    @assert length(_partition[1]) == deriv_num + 1 "partition should have $deriv_num+1 entries"

    diagpara = Vector{DiagPara}()
    MaxOrder = maximum([p[1] for p in _partition])
    MinOrder = minimum([p[1] for p in _partition])
    dict_graphs = Dict{NTuple{deriv_num + 1,Int},Vector{Graph}}()
    spinPolarPara = 2.0 / paramc.spin - 1

    for order in MinOrder:MaxOrder
        partition_ind = findall(p -> p[1] == order, _partition)
        partition_order = [_partition[i] for i in partition_ind]

        diagrams = GV.diagsGV(diagtype, order, spinPolarPara=spinPolarPara, filter=filter)
        optimize!(diagrams, level=optimize_level)
        optimize!(diagrams, level=optimize_level)

        renormalization_orders = Int[]
        for i in 1:deriv_num
            push!(renormalization_orders, maximum([p[i+1] for p in partition_order]))
        end
        # renormalization_orders = [Max_GD_o, Max_ID_o, extra_deriv_orders...]
        dict_graph_order = taylorAD(diagrams, renormalization_orders, leaf_dep_funcs)
        for key in keys(dict_graph_order)
            p = (order, key...)
            if p in _partition
                dict_graphs[p] = dict_graph_order[key]
            end
        end
    end

    extT_labels = Vector{Vector{Int}}[]
    partition = sort(collect(keys(dict_graphs)))
    for p in partition
        push!(diagpara, diagPara(_diagtype(diagtype), paramc.isDynamic, p[1], paramc.spin, filter))
        push!(extT_labels, [collect(g.properties.extT) for g in dict_graphs[p]])
    end
    return (partition, diagpara, dict_graphs, extT_labels)
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