using CodecZlib
using ElectronLiquid
using JLD2

"""
    function remove_para_from_jld2_data(filename; save=false, log=true)

Removes redundant entries of type ParaMC from a JLD2 archive.
For example, a data entry `(a, p::ParaMC, b, c)` with key `short(p)` will be replaced by `(a, b, c)`.

# Arguments
- `filename`: the name of the JLD2 archive to be updated
- `save`: whether to write the updated archive to file `filename`
- `log`: whether write log info to file `remove_para_from_jld2_data.log`
"""
function remove_para_from_jld2_data(filename; save=false, log=false)
    io = IOBuffer()
    println.([io, stdout], "\nChecking for redundant ParaMC in $filename...\n")

    has_updates = false
    data = load(filename)
    for k in keys(data)
        # Convert JLD2 entries like (a, p::ParaMC, b, c, ...) to (a, b, c, ...)
        if data[k] isa Tuple
            idx_para = findall(v -> v isa ParaMC, data[k])
            if isempty(idx_para) == false
                has_updates = true
                if save
                    println.([io, stdout], "Removing ParaMC in data with key $k...\n")
                    data[k] = Tuple(v for (i, v) in enumerate(data[k]) if i ∉ idx_para)
                else
                    println.([io, stdout], "Will remove ParaMC in data with key $k...\n")
                end
            else
                println.([io, stdout], "No ParaMC in data with key $k...\n")
            end
        end
    end

    if has_updates == false
        println.([io, stdout], "No changes needed to $filename.\n")
        return
    end

    if save
        # Backup the JLD2 archive
        suffix = 0
        backup_name = "$(filename).bak"
        if isfile(backup_name)
            while isfile(backup_name)
                suffix += 1
                backup_name = "$(filename).bak$(suffix)"
            end
        end
        println.([io, stdout], "Creating a backup at $(backup_name)...\n")
        cp(filename, backup_name)

        # Save the updated JLD2 archive
        println.([io, stdout], "Saving data...\n")
        jldopen(filename, "w"; compress=true) do f
            for (k, v) in data
                f[k] = v
            end
        end
    end

    if log
        # Save the log info to file
        open("remove_para_from_jld2_data.log", "a+") do f
            write(f, String(take!(io)))
        end
    end

    return
end
