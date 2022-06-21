using JLD2
using Measurements
# using Plots
include("../common/parameter.jl")
include("../common/counterterm.jl")
using .CounterTerm


function zfactor(idata)
    return (idata[2] - idata[1]) / (2π / β)
end

function loaddata(FileName)
    f = jldopen(FileName, "r")
    avg, std = f["avg"], f["std"]
    order = f["order"]
    _partition = f["partition"]
    rdata, idata = Dict(), Dict()
    for (ip, p) in enumerate(_partition)
        rdata[p] = measurement.(real(avg[ip, :, :]), real(std[ip, :, :]))
        idata[p] = measurement.(imag(avg[ip, :, :]), imag(std[ip, :, :]))
    end
    return order, _partition, rdata, idata
end

if abspath(PROGRAM_FILE) == @__FILE__

    @assert length(ARGS) >= 1 "One argument for the data file name is required!"
    filename = ARGS[1]
    _order, _partition, rdata, idata = loaddata(filename)

    rdata = mergeInteraction(rdata)
    idata = mergeInteraction(idata)

    df = fromFile()
    δμ = muCT(df, paraid, _order)
    println(δμ)

end