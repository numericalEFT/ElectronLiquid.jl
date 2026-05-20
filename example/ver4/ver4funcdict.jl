using FeynmanDiagram, ElectronLiquid

const isDynamic = false

function func_dict_str(o)
    str = ""
    str *= "const evalfuncParquetAD_map = Dict(\n"
    # for order in 1:o
    order = o
    para = UEG.ParaMC(rs=1.0, beta=25, order=order, isDynamic=isDynamic)
    _partition = UEG.partition(order, offset=0)
    partition = Vector{Tuple{Int64,Int64,Int64}}()
    for (o, sOrder, vOrder) in _partition
        o == 0 && sOrder > 0 && continue
        push!(partition, (o, sOrder, vOrder))
    end
    if isDynamic
        diagram = Ver4.diagramParquetDynamic_load(para, partition)
    else
        diagram = Ver4.diagramParquet_load(para, partition)
    end
    _partition = diagram[1]
    println(partition)
    println(_partition)
    for p in _partition
        str *= "($order, $(p[1]), $(p[2]), $(p[3])) => eval_ver4O$(order)ParquetAD$(p[1])$(p[2])$(p[3])!, \n"
    end
    # end
    str *= ")\n"
    return str
end

fname = "../../src/vertex4/source_codeParquetAD/func_dict_ParquetAD.jl"
# fname = "./src/vertex4/source_codeParquetAD/dynamic/func_dict_ParquetADDynamic.jl"

open(fname, "w") do f
    write(f, func_dict_str(1))
end
