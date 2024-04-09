function compileC_ParquetAD_toFiles(FeynGraphs, maxMomNum::Int, diagname::String;
    datatype::DataType=Float64, root_dir=joinpath(@__DIR__, "source_codeParquetAD"),
    c_source=joinpath(root_dir, "func_$(diagname)_ParquetAD.c"),
    lib_path=root_dir, lib_name="$(diagname)_ParquetAD", compiler::String="gcc", isnative::Bool=false
)
    noresponse = true
    has_extT = true
    if length(FeynGraphs) == 3
        partition, diagpara, FeynGraphs = FeynGraphs
        has_extT = false
    elseif length(FeynGraphs) == 4
        partition, diagpara, FeynGraphs, extT_labels = FeynGraphs
    elseif length(FeynGraphs) == 5
        partition, diagpara, FeynGraphs, extT_labels, spin_conventions = FeynGraphs
        noresponse = false
    else
        error("Invalid FeynGraphs format")
    end

    ### compile the Parquet + Taylor-AD generated Graphs to C language source code
    println("compiling to c code")
    leaf_maps = ParquetADcompileC_toFile(partition, FeynGraphs, c_source, datatype=datatype)
    println("compiled")

    # ### compile the C language 
    println("compiling c code to .so lib")
    ParquetADcompileC_so(partition, c_source, lib_path, lib_name, datatype=datatype, compiler=compiler, isnative=isnative)
    println("compiled")
    ParquetAD_Cwrapper(partition, lib_path, lib_name, datatype=datatype)

    println("saving other info")
    ### save the leafs information and the loopbasis to CSV files
    leafinfo_toFile(partition, leaf_maps, maxMomNum, root_dir, diagname)

    ### save the external tau variables' indexes and spin channel to a jld2 file

    if has_extT
        if noresponse
            extvar_toFile(partition, root_dir, diagname, extT_labels)
        else
            extvar_toFile(partition, root_dir, diagname, extT_labels, spin_conventions)
        end
    end
    println("saved")
end

function ParquetADcompileC_toFile(partition, FeynGraphs, c_source::String; datatype::DataType=Float64)

    ### compile and save the generated Cfunction to a source code file.
    # leaf_maps = Vector{Dict{Int,Union{Graph,FeynmanGraph}}}()
    leaf_maps = Vector{Dict{Int,Graph}}()
    for key in partition
        key_str = join(string.(key))
        leafmap = Compilers.compile_C(FeynGraphs[key], c_source, func_name="eval_graph$(key_str)", datatype=datatype)
        push!(leaf_maps, leafmap)
    end

    return leaf_maps
    # ### save the leafs information to CSV files 
    # leafinfo_toFile(partition, leaf_maps, labelProd, root_dir)
end

function extvar_toFile(partition, root_dir::String, diagname::String, vars...)
    jldopen(joinpath(root_dir, "extvars_$(diagname).jld2"), "w") do f
        for (i, key) in enumerate(partition)
            key_str = join(string.(key))
            f[key_str] = Tuple([v[i] for v in vars])
        end
    end
end

function ParquetADcompileC_so(partition, c_source::String, lib_path::String, lib_name::String;
    datatype::DataType=Float64, compiler::String="gcc", isnative::Bool=false
)
    lib = isnative ? joinpath(lib_path, "$(lib_name)_native.so") : joinpath(lib_path, "$lib_name.so")
    ### compile C source file to *.so library
    if isnative
        command = `$compiler -shared -fPIC -o $lib $c_source -O3 -march=native`
    else
        command = `$compiler -shared -fPIC -o $lib $c_source -O3`
    end
    run(command)

    ### generate the C wrapper Julia functions and save them in the same path of C library
    ParquetAD_Cwrapper(partition, lib_path, lib_name, datatype=datatype)
end

function ParquetAD_Cwrapper(partition, lib_path::String, lib_name::String; datatype::DataType=Float64)
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

function leafinfo_toFile(partition, leaf_maps::Vector{Dict{Int,Graph}}, maxMomNum::Int,
    root_dir::String, diagname::String)
    leafStates, loopbasis = FeynmanDiagram.leafstates(leaf_maps, maxMomNum)
    len = length(leafStates)

    for (ikey, key) in enumerate(partition)
        key_str = join(string.(key))
        df = DataFrame([leafStates[idx][ikey] for idx in 1:len], :auto)
        CSV.write(joinpath(root_dir, "leafinfo_$(diagname)_$key_str.csv"), df)
    end

    ### save the loop basis to a CSV file for the maximum order
    order = maximum(p[1] for p in partition)
    df = DataFrame(loopbasis, :auto)
    CSV.write(joinpath(root_dir, "loopBasis_$(diagname)_maxOrder$order.csv"), df)
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