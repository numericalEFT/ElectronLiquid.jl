function compileC_ParquetAD_toFiles(order, partition, FeynGraphs, maxloopNum::Int; datatype::DataType=Float64,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD"), c_source=joinpath(root_dir, "func_O$(order)_ver4ParquetAD.c"),
    lib_path=root_dir, lib_name="ver4O$(order)ParquetAD", compiler::String="gcc", isnative::Bool=false)

    ### compile the Parquet + Taylor-AD generated Graphs to C language source code
    println("compiling to c code")
    leaf_maps = ParquetADcompileC_toFile(order, partition, FeynGraphs, root_dir, c_source=c_source)
    println("compiled")

    # ### compile the C language 
    # println("compiling c code to .so lib")
    # ParquetADcompileC_so(order, partition, datatype; c_source=c_source,
    #     lib_path=lib_path, lib_name=lib_name, compiler=compiler, isnative=isnative)
    # println("compiled")
    ParquetAD_Cwrapper(order, partition, datatype, lib_path=lib_path, lib_name=lib_name)

    println("saving other info")
    ### save the leafs information and the loopbasis to CSV files
    leafinfo_toFile(order, partition, leaf_maps, maxloopNum, root_dir)

    ### save the external tau variables' indexes and spin channel to a jld2 file
    extT_and_spin_toFile(order, partition, FeynGraphs, root_dir)
    println("saved")
end

function ParquetADcompileC_toFile(order, partition, FeynGraphs,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD"); c_source=joinpath(root_dir, "func_O$(order)_ver4ParquetAD.c"), datatype::DataType=Float64)

    ### compile and save the generated Cfunction to a source code file.
    # leaf_maps = Vector{Dict{Int,Union{Graph,FeynmanGraph}}}()
    leaf_maps = Vector{Dict{Int,Graph}}()
    for key in partition
        key_str = join(string.(key))
        leafmap = Compilers.compile_C(FeynGraphs[key][1], c_source; func_name="eval_graph$(key_str)", datatype=datatype)
        push!(leaf_maps, leafmap)
    end

    return leaf_maps
    # ### save the leafs information to CSV files 
    # leafinfo_toFile(partition, leaf_maps, labelProd, root_dir)
end

function extT_and_spin_toFile(order, partition, FeynGraphs, root_dir=joinpath(@__DIR__, "source_codeParquetAD"))
    jldopen(joinpath(root_dir, "extT_spin_O$(order)_ParquetAD.jld2"), "w") do f
        for key in partition
            key_str = join(string.(key))
            f[key_str] = (FeynGraphs[key][2], FeynGraphs[key][3])
        end
    end
end

function ParquetADcompileC_so(order, partition, datatype::DataType=Float64; c_source=joinpath(@__DIR__, "source_codeGV", "func_O$(order)_ver4ParquetAD.c"),
    lib_path=joinpath(@__DIR__, "source_codeParquetAD"), lib_name="ver4O$(order)ParquetAD", compiler::String="gcc", isnative::Bool=false)

    lib = isnative ? joinpath(lib_path, "$(lib_name)_native.so") : joinpath(lib_path, "$lib_name.so")
    ### compile C source file to *.so library
    if isnative
        command = `$compiler -shared -fPIC -o $lib $c_source -O3 -march=native`
    else
        command = `$compiler -shared -fPIC -o $lib $c_source -O3`
    end
    run(command)

    ### generate the C wrapper Julia functions and save them in the same path of C library
    ParquetAD_Cwrapper(order, partition, datatype, lib_path=lib_path, lib_name=lib_name)
end

function ParquetAD_Cwrapper(order, partition, datatype::DataType=Float64; lib_path=joinpath(@__DIR__, "source_codeParquetAD"), lib_name="ver4O$(order)ParquetAD")
    ctype_str = julia_to_C_aliasstr(datatype)
    str = ""
    for key in partition
        key_str = join(string.(key))
        str *= "\nfunction eval_$(lib_name)$(key_str)!(root::Vector{$datatype}, leafVal::Vector{$datatype})\n"
        str *= "    @ccall joinpath(@__DIR__, \"$lib_name.so\").eval_graph$(key_str)(root::Ptr{$ctype_str}, leafVal::Ptr{$ctype_str})::Cvoid\n"
        str *= "end"
    end

    ### save the C wrapper Julia functions in the same path of C library
    open(joinpath(lib_path, "Cwrapper_$lib_name.jl"), "w") do f
        write(f, str)
    end
end

function leafinfo_toFile(order, partition, leaf_maps::Vector{Dict{Int,Graph}}, maxloopNum::Int, root_dir=joinpath(@__DIR__, "source_codeParquetAD/"))
    leafStates, loopbasis = FeynmanDiagram.leafstates_diagtree(leaf_maps, maxloopNum)
    len = length(leafStates)

    for (ikey, key) in enumerate(partition)
        key_str = join(string.(key))
        df = DataFrame([leafStates[idx][ikey] for idx in 1:len], :auto)
        CSV.write(joinpath(root_dir, "leafinfo_O$(order)_Parquet$key_str.csv"), df)
    end

    ### save the loop basis to a CSV file for the maximum order
    df = DataFrame(loopbasis, :auto)
    CSV.write(joinpath(root_dir, "loopBasis_ParquetADmaxOrder$order.csv"), df)
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