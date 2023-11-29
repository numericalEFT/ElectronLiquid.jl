function GVcompileJulia_toFile(maxOrder::Int, FeynGraphs, labelProd::LabelProduct, root_dir=joinpath(@__DIR__, "source_codeGV"))
    leaf_maps = Vector{Dict{Int,FeynmanGraph}}()
    ### compile and save the generated Julia function to a source code file.
    for order in 1:maxOrder
        partition = UEG.partition_order(order)
        for key in partition
            key_str = join(string.(key))
            leafmap = Compilers.compile_Julia(FeynGraphs[key][1], joinpath(root_dir, "func_sigmaGV_o$order.jl"); func_name="eval_graph$(key_str)!")
            push!(leaf_maps, leafmap)
        end
    end

    ### save the leafs information to a CSV file
    leafinfo_toFile(maxOrder, leaf_maps, labelProd, root_dir)
end

function GVcompileC_toFile(partition, FeynGraphs, labelProd::LabelProduct, 
    root_dir=joinpath(@__DIR__, "source_codeGV"); datatype::DataType=Float64)

    ### compile and save the generated Cfunction to a source code file.
    leaf_maps = Vector{Dict{Int,FeynmanGraph}}()
    for key in partition
        key_str = join(string.(key))
        leafmap = Compilers.compile_C(FeynGraphs[key][1], joinpath(root_dir, "func_sigmaGV.c"); func_name="eval_graph$(key_str)", datatype=datatype)
        push!(leaf_maps, leafmap)
    end

    ### save the leafs information to CSV files 
    leafinfo_toFile(partition, leaf_maps, labelProd, root_dir)
end

function leafinfo_toFile(partition, leaf_maps::Vector{Dict{Int,FeynmanGraph}}, labelProd::LabelProduct, root_dir=joinpath(@__DIR__, "source_codeGV"))
    leafStates = FeynmanDiagram.leafstates(leaf_maps, labelProd)
    len = length(leafStates)
    for (ikey, key) in enumerate(partition)
        key_str = join(string.(key))
        df = DataFrame([leafStates[idx][ikey] for idx in 1:len], :auto)
        CSV.write(joinpath(root_dir, "leafinfo_GV$key_str.csv"), df)
    end
end

function loopbasis_toFile(maxOrder::Int, labelProd::LabelProduct, root_dir=joinpath(@__DIR__, "source_codeGV"))
    ### save the loop basis to a CSV file for the given maximum order
    df = DataFrame(labelProd.labels[end], :auto)
    CSV.write(joinpath(root_dir, "loopBasis_GVmaxOrder$maxOrder.csv"), df)
end

function GVcompileC_so(partition, datatype::DataType=Float64; c_source=joinpath(@__DIR__, "source_codeGV", "func_sigmaGV.c"), 
    lib_path=joinpath(@__DIR__, "source_codeGV"), lib_name="sigmaGV", compiler::String="gcc")
    ### compile C source file to *.so library
    lib = joinpath(lib_path, "$lib_name.so")
    run(`$compiler -shared -fPIC -o $lib $c_source -O3`)

    ### generate the C wrapper Julia functions and save them in the same path of C library
    GV_Cwrapper(partition, datatype, lib_path=lib_path, lib_name=lib_name)
end

function GV_Cwrapper(partition, datatype::DataType=Float64; lib_path=joinpath(@__DIR__, "source_codeGV"), lib_name="sigmaGV")
    ctype_str = julia_to_C_aliasstr(datatype)
    str=""
    for key in partition
        key_str = join(string.(key))
        str *= "\nfunction eval_graph$(key_str)!(root::Vector{$datatype}, leafVal::Vector{$datatype})\n"
        str *= "    @ccall joinpath(@__DIR__, \"$lib_name.so\").eval_graph$(key_str)(root::Ptr{$ctype_str}, leafVal::Ptr{$ctype_str})::Cvoid\n"
        str *= "end"
    end

    ### save the C wrapper Julia functions in the same path of C library
    open(joinpath(lib_path, "Cwrapper_$lib_name.jl"), "w") do f 
        write(f, str)
    end
end

function julia_to_C_aliasstr(type::DataType)
    if type == Float64
        return "Cdouble"
    elseif type == Float32
        return "Cfloat"
    elseif type == Int64
        return "Clonglong"
    elseif type == Int32
        return "Cint"
    elseif type == ComplexF32
        return "ComplexF32"
    elseif type == ComplexF64
        return "ComplexF64"
    elseif type <: Array
        return julia_to_C_aliasstr(eltype(type)) * "*"
    else
        error("Unsupported type")
    end
end