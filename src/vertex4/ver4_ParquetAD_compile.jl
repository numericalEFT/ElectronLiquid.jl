function compileC_ParquetAD_toFiles(partition, FeynGraphs, maxloopNum::Int; datatype::DataType=Float64,
    root_dir=joinpath(@__DIR__, "source_codeGV"), c_source=joinpath(root_dir, "func_sigmaParquetAD.c"),
    lib_path=root_dir, lib_name="sigmaParquetAD", compiler::String="gcc", isnative::Bool=false)

    ### compile the Parquet + Taylor-AD generated Graphs to C language source code
    leaf_maps = GVcompileC_toFile(partition, FeynGraphs, root_dir, c_source=c_source)

    ### compile the C language 
    GVcompileC_so(partition, datatype; c_source=c_source,
        lib_path=lib_path, lib_name=lib_name, compiler=compiler, isnative=isnative)

    ### save the leafs information and the loopbasis to CSV files
    leafinfo_toFile(partition, leaf_maps, maxloopNum, root_dir)

    ### save the external tau variables' indexes to a jld2 file
    extT_toFile(partition, FeynGraphs, root_dir)
end

function GVcompileC_toFile(partition, FeynGraphs,
    root_dir=joinpath(@__DIR__, "source_codeGV"); c_source=joinpath(root_dir, "func_sigmaGV.c"), datatype::DataType=Float64)

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

function extT_toFile(partition, FeynGraphs, root_dir=joinpath(@__DIR__, "source_codeGV"))
    jldopen(joinpath(root_dir, "extT_ParquetAD.jld2"), "w") do f
        for key in partition
            key_str = join(string.(key))
            f[key_str] = FeynGraphs[key][2]
        end
    end
end

function GVcompileC_so(partition, datatype::DataType=Float64; c_source=joinpath(@__DIR__, "source_codeGV", "func_sigmaGV.c"),
    lib_path=joinpath(@__DIR__, "source_codeGV"), lib_name="sigmaGV", compiler::String="gcc", isnative::Bool=false)

    lib = isnative ? joinpath(lib_path, "$(lib_name)_native.so") : joinpath(lib_path, "$lib_name.so")
    ### compile C source file to *.so library
    if isnative
        command = `$compiler -shared -fPIC -o $lib $c_source -O3 -march=native`
    else
        command = `$compiler -shared -fPIC -o $lib $c_source -O3`
    end
    run(command)

    ### generate the C wrapper Julia functions and save them in the same path of C library
    GV_Cwrapper(partition, datatype, lib_path=lib_path, lib_name=lib_name)
end

function GV_Cwrapper(partition, datatype::DataType=Float64; lib_path=joinpath(@__DIR__, "source_codeGV"), lib_name="sigmaGV")
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