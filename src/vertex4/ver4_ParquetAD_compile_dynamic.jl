
# compile ver4 diagrams with dynamic interaction
# only difference is storing more leafinfo

function compileC_ParquetAD_toFiles_dynamic(order, partition, FeynGraphs, maxloopNum::Int; datatype::DataType=Float64,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/dynamic"), c_source=joinpath(root_dir, "func_O$(order)_ver4ParquetADDynamic.c"),
    lib_path=root_dir, lib_name="ver4O$(order)ParquetADDynamic", compiler::String="gcc", isnative::Bool=false)

    ### compile the Parquet + Taylor-AD generated Graphs to C language source code
    println("compiling to c code")
    leaf_maps = ParquetADcompileC_toFile(order, partition, FeynGraphs, root_dir, c_source=c_source)
    println("compiled")

    println("saving other info")
    ### save the leafs information and the loopbasis to CSV files
    leafinfo_toFile_dynamic(order, partition, leaf_maps, maxloopNum, root_dir)

    # ### compile the C language 
    println("compiling c code to .so lib")
    ParquetADcompileC_so(order, partition, datatype; c_source=c_source,
        lib_path=lib_path, lib_name=lib_name, compiler=compiler, isnative=isnative)
    println("compiled")
    ParquetAD_Cwrapper(order, partition, datatype, lib_path=lib_path, lib_name=lib_name)

    ### save the external tau variables' indexes and spin channel to a jld2 file
    extT_and_spin_toFile(order, partition, FeynGraphs, root_dir)
    println("saved")
end

function leafinfo_toFile_dynamic(order, partition, leaf_maps::Vector{Dict{Int,Graph}}, maxloopNum::Int, root_dir=joinpath(@__DIR__, "source_codeParquetAD/dynamic/"))
    leafStates, loopbasis = FeynmanDiagram.leafstates(leaf_maps, maxloopNum)
    len = length(leafStates)

    leafstates = Vector{Vector{LeafStateADDynamic}}()
    leafvalues = Vector{Vector{Float64}}()

    for (ikey, key) in enumerate(partition)
        key_str = join(string.(key))
        df = DataFrame([leafStates[idx][ikey] for idx in 1:len], :auto)
        leafstates_par = Vector{LeafStateADDynamic}()
        # for row in eachrow(df)
        println(size(df), key)
        for idx in 1:size(df)[1]
            row = df[idx, :]
            # row = df[:, idx]
            diagid = leaf_maps[ikey][idx].properties
            tau_num = interactionTauNum(diagid.para)
            idorder = row[3]
            if row[2] == 2
                # interaction line
                if diagid.type == Instant
                    row[2] = 2
                elseif diagid.type == Dynamic
                    row[2] = 4
                end
                idorder = diagid.order
            end
            # println(row)
            # @assert diagid.order[1:2] == row[3] "$(diagid.order[1:2]) == $(row[3])"
            push!(leafstates_par, LeafStateADDynamic(row[2], idorder, row[4:end]..., tau_num))
        end
        push!(leafstates, leafstates_par)
        push!(leafvalues, df[!, names(df)[1]])
    end

    jldopen(joinpath(root_dir, "leafinfo_O$(order).jld2"), "w") do f
        f["leafstates"] = leafstates
        f["values"] = leafvalues
        println(leafvalues)
    end

    ### save the loop basis to a CSV file for the maximum order
    df = DataFrame(loopbasis, :auto)
    CSV.write(joinpath(root_dir, "loopBasis_ParquetADmaxOrder$order.csv"), df)
end
