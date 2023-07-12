using ElectronLiquid
using Measurements
using JLD2

const filename = "data_K.jld2"

function mu(data)
    return real(data[1, 1])
end

function process(FileName, isSave)
    f = jldopen(FileName, "r")
    for key in keys(f)
        println(key)
        value = f[key]
        para = UEG.ParaMC(key)
        ngrid, kgrid, data = value[1], value[2], value[3]
        dim, β, kF, me = para.dim, para.β, para.kF, para.me
        kF_label = searchsortedfirst(kgrid, kF)
        println(kF_label)
        kF_half_label = searchsortedfirst(kgrid, kF / 2)
        for (p, val) in data
            ek = real.(val[1, :] .- val[1, 1]) ./ (kgrid .^ 2 / (2me))
            println("order: $p, mass1 = $(ek[kF_label]), mass2 = $(ek[kF_half_label])")
            println(ek)
        end

    end
end

if abspath(PROGRAM_FILE) == @__FILE__

    # @assert length(ARGS) >= 1 "One argument for the data file name is required!"
    # filename = ARGS[1]
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end
    process(filename, isSave)
end