import ElectronLiquid: UEG, CounterTerm
using ElectronGas
using Measurements
using Printf
using JLD2
import FeynmanDiagram: SpinSpin, ChargeCharge

dim = 3
rs = [4.0,]
mass2 = [1.0,]
Fs = [-0.0,]
beta = [25.0,]
order = [4,]
isDynamic = false
# response = ChargeCharge
response = SpinSpin

const filenameN = "data_n.jld2"
const filenameP = "data_$(response)Polar_K.jld2"

# const parafilename = "para.csv"
# const parafilename = "para_wn_1minus0.csv"

function process(para, datatuple, isSave)
    # df = CounterTerm.fromFile(parafilename)
    data = datatuple[1]
    printstyled(UEG.short(para), color=:yellow)
    println()
    # for p in sort([k for k in keys(data)])
    #     println("$p: μ = $(data[p])")
    # end
    _mu = Dict()
    for (p, val) in data
        _mu[p] = val
    end

    # dzi, dmu, dz = CounterTerm.sigmaCT(para.order, _mu, _z)
    dmu = CounterTerm.densityCT(para.order, _mu, verbose=1)
    for i in eachindex(dmu)
        println("dmu[$i]: ", dmu[i])
    end

    return dmu
end

function loaddata(para, FileName)
    key = UEG.short(para)
    f = jldopen(FileName, "r")

    p = UEG.ParaMC(key)
    ngrid, kgrid, data = f[key]
    # println(data)
    _partition = UEG.partition(para.order)
    partition = Vector{eltype(_partition)}()
    for (order, sOrder, vOrder) in _partition
        order == 1 && vOrder > 0 && continue
        push!(partition, (order, sOrder, vOrder))
    end
    rdata, idata = Dict(), Dict()
    for (ip, p) in enumerate(partition)
        rdata[p] = real(data[p][:, :])
        idata[p] = imag(data[p][:, :])
    end
    return ngrid, kgrid, rdata, idata
end

function derive_polar(para, dmu)
    ngrid, kgrid, rdata, idata = loaddata(para, filenameP)
    # println("Re part:")
    # for (p, val) in rdata
    #     println(p, val)
    # end
    # println("Im part:")
    # for (p, val) in idata
    #     println(p, val)
    # end
    rPolar = CounterTerm.chemicalpotential_renormalization(para.order, rdata, dmu)
    iPolar = CounterTerm.chemicalpotential_renormalization(para.order, idata, dmu)
    return ngrid, kgrid, rPolar, iPolar
end

function derive_Kxc(para, ngrid, kgrid, rPolar, iPolar)
    Kxc = similar(rPolar[1])
    Gfactor = similar(rPolar[1])
    Π0 = zeros(length(ngrid), length(kgrid))
    for (i, n) in enumerate(ngrid)
        for (j, q) in enumerate(kgrid)
            Π0[i, j] = ElectronGas.Polarization.Polarization0_FiniteTemp(q, n, para.basic) * para.spin * para.massratio
        end
    end
    Π0 = measurement.(Π0, 0.0)

    println(rPolar)
    rPolar[1] = -Π0

    order = length(rPolar)
    for o in 2:order
        Π = sum(-rPolar[1:o])
        println("Kxc Order $(o-1): ", 1 ./ Π0 - 1 ./ Π)
    end

    Π = sum(-rPolar)
    Kxc = 1 ./ Π0 - 1 ./ Π
    println("Kxc:", Kxc)
    for i in eachindex(ngrid)
        Gfactor[i, :] = -Kxc[i, :] .* kgrid .* kgrid / 8π
    end
    println("G:", Gfactor)
    return Kxc, Gfactor
end

# function Kxc_renormalization(order, data)

# end

function derive_Kxc_renorm(para, ngrid, kgrid, rPolar, iPolar)
    println()
    Π0 = zeros(length(ngrid), length(kgrid))
    for (i, n) in enumerate(ngrid)
        for (j, q) in enumerate(kgrid)
            Π0[i, j] = ElectronGas.Polarization.Polarization0_FiniteTemp(q, n, para.basic) * para.spin * para.massratio
        end
    end
    Π0 = measurement.(Π0, 0.0)
    println("Π0:", Π0)

    rPolar[1] = -Π0
    rPolar = -rPolar
    for o in eachindex(rPolar)
        println("Π Order $o: ", rPolar[o])
    end

    order = length(rPolar) - 1
    Kxc = [similar(Π0) for _ in 1:order]
    if order >= 1
        Kxc[1] = rPolar[2] ./ Π0 .^ 2
    end
    if order >= 2
        Kxc[2] = rPolar[3] ./ Π0 .^ 2 .- rPolar[2] .^ 2 ./ Π0 .^ 3
    end
    if order >= 3
        Kxc[3] = rPolar[4] ./ Π0 .^ 2 - 2rPolar[2] .* rPolar[3] ./ Π0 .^ 3 + rPolar[2] .^ 3 ./ Π0 .^ 4
    end
    if order >= 4
        Kxc[4] = rPolar[5] ./ Π0 .^ 2 - 2rPolar[2] .* rPolar[4] ./ Π0 .^ 3 - rPolar[3] .^ 2 ./ Π0 .^ 3
        +3rPolar[3] .* rPolar[2] .^ 2 ./ Π0 .^ 4 - rPolar[2] .^ 4 ./ Π0 .^ 5
    end
    for o in 1:order
        println("Kxc Order $o: ", Kxc[o])
    end
    # println(Kxc[1] .* Π0 .^ 2)
    # println(Kxc[2] .* Π0 .^ 2 .+ Kxc[1] .^ 2 .* Π0 .^ 3)
    # println(Kxc[3] .* Π0 .^ 2 .+ 2Kxc[1] .* Kxc[2] .* Π0 .^ 3 .+ Kxc[1] .^ 3 .* Π0 .^ 4)

    Kxc_sum = sum(Kxc)
    println("sum Kxc:", Kxc_sum)
    Gfactor = similar(Π0)
    for i in eachindex(ngrid)
        Gfactor[i, :] = -Kxc_sum[i, :] .* kgrid .* kgrid / 8π
    end
    println("G:", Gfactor)
    return Kxc, Gfactor
end

if abspath(PROGRAM_FILE) == @__FILE__

    # @assert length(ARGS) >= 1 "One argument for the data file name is required!"
    # filename = ARGS[1]
    isSave = false
    if length(ARGS) >= 1 && (ARGS[1] == "s" || ARGS[1] == "-s" || ARGS[1] == "--save" || ARGS[1] == " save")
        # the second parameter may be set to save the derived parameters
        isSave = true
    end

    f_n = jldopen(filenameN, "r")

    for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
        para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim)
        kF = para.kF
        for key in keys(f_n)
            loadpara = UEG.ParaMC(key)
            if UEG.paraid(loadpara) == UEG.paraid(para)
                dmu = process(para, f_n[key], isSave)
                ngrid, kgrid, rPolar, iPolar = derive_polar(para, dmu)
                derive_Kxc(para, ngrid, kgrid, rPolar, iPolar)
                derive_Kxc_renorm(para, ngrid, kgrid, rPolar, iPolar)
            end
        end
    end
end
