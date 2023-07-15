using CodecZlib
using ElectronLiquid
using JLD2

# Converts JLD2 data from old to new ParaMC format on load by adding the `initialized` field (see: https://juliaio.github.io/JLD2.jl/stable/advanced/)
# NOTE: Requires the type name `ElectronLiquid.UEG.ParaMC` to be explicitly specified
JLD2.rconvert(::Type{ElectronLiquid.UEG.ParaMC}, nt::NamedTuple) = ElectronLiquid.UEG.ParaMC(; nt..., initialized=false)

"""
    function update_jld2_data(filename; save=false, log=true)

Removes redundant entries of type ParaMC from a JLD2 archive and upgrades key formats (field `order` now present in paraid).
For example, the data entry `(a, p::ParaMC, b, c)` will be replaced by `(a, b, c)`, and its key will be updated to include
the substring `order_\$(p.order)`.

# Arguments
- `filename`: the name of the JLD2 archive to be updated
- `save`: whether to write the updated archive to file `filename`
- `log`: whether to write log info to file `update_jld2_data.log`
"""
function update_jld2_data(filename; save=false, log=false)
    io = IOBuffer()
    println.([io, stdout], "\nChecking for redundant ParaMC entries in $filename...\n")

    # Load the JLD2 archive. 
    # NOTE: A typemap must be specified to update outdated ParaMC types stored in the data to match the new struct definition.
    data = load(filename; typemap=Dict("ElectronLiquid.UEG.ParaMC" => JLD2.Upgrade(ElectronLiquid.UEG.ParaMC)))
    updated_data = Dict()

    # Convert JLD2 data entries like (a, p::ParaMC, b, c, ...) to (a, b, c, ...)
    has_updates = false
    for k in keys(data)
        new_key = k
        try
            # First, check that we can parse this key as a ParaMC object
            UEG.ParaMC(k)
            # Add field `order` to keys and remove all ParaMC instances
            if data[k] isa Tuple
                idx_para = findfirst(v -> v isa ParaMC, data[k])
                if isnothing(idx_para)
                    println.([io, stdout], "No ParaMC in data entry with key $k...\n")
                else
                    has_updates = true
                    # Get updated key name including field `order`
                    new_key = UEG.short(data[k][idx_para])
                    if save
                        println.([io, stdout], "Removing ParaMC in data entry with key $k...")
                        println.([io, stdout], "Updating key name to $new_key...\n")
                        # Update this entry
                        data[k] = Tuple(v for (i, v) in enumerate(data[k]) if i != idx_para)
                    else
                        println.([io, stdout], "Will remove ParaMC in data entry with key $k...")
                        println.([io, stdout], "Will update key name to $new_key...\n")
                    end
                end
            end
        catch
            println.([io, stdout], "Failed to parse key $k as a ParaMC object, skipping it...\n")
        end
        updated_data[new_key] = data[k]
    end
    if has_updates == false
        println.([io, stdout], "No changes needed to $filename.\n")
        println.([io, stdout], "Done!\n")
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
            for (k, v) in updated_data
                f[k] = v
            end
        end
    end
    println.([io, stdout], "Done!\n")

    if log
        # Save the log info to file
        open("update_jld2_data.log", "a+") do f
            write(f, String(take!(io)))
        end
    end
    return
end
