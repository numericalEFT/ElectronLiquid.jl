module UEG
using StaticArrays
using Parameters
using ElectronGas: Parameter
using CompositeGrids

export INL, INR, OUTL, OUTR, DI, EX
export ParaMC, Weight

@with_kw struct ParaMC
    ### fundamental parameters
    beta::Float64
    rs::Float64
    order::Int = 2
    Fs::Float64 = -0.0
    Fa::Float64 = -0.0

    mass2::Float64 = 1e-5
    massratio::Float64 = 1.0
    dim::Int = 3
    spin::Int = 2

    #### predefined #################################
    isFock = false
    isDynamic = true # if true, use dynamic interaction, otherwise use instant interaction

    #### derived parameters ###########
    basic = Parameter.rydbergUnit(1.0 / beta, rs, dim, Λs=mass2)

    kF = basic.kF
    EF = basic.EF
    β = basic.β
    maxK = 6 * basic.kF
    me = basic.me
    ϵ0 = basic.ϵ0
    e0 = basic.e0

    fs = Fs / basic.NF / massratio
    fa = Fa / basic.NF / massratio

    ##########   effective interaction and counterterm ###############
    qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, maxK], [0.0, 2kF], 16, 0.01 * kF, 8)
    τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)

    dW0 = KO(basic, qgrid, τgrid, mass2, massratio, Fs, Fa)
    cRs = [counterKO(basic, qgrid, τgrid, o, mass2, massratio, Fs, Fa) for o in 1:order]
end

paraid(p::ParaMC) = Dict(
    "dim" => p.dim,
    "rs" => p.rs,
    "beta" => p.beta,
    "mass2" => p.mass2,
    "Fs" => p.Fs,
    "Fa" => p.Fa,
    "massratio" => p.massratio,
    "spin" => p.spin,
    "isFock" => p.isFock,
    "isDynamic" => p.isDynamic
)

const INL, OUTL, INR, OUTR = 1, 2, 3, 4
const DI, EX = 1, 2

mutable struct Weight <: FieldVector{2,Float64}
    d::Float64
    e::Float64
    Weight() = new(0.0, 0.0)
    Weight(d, e) = new(d, e)
end

const Base.zero(::Type{Weight}) = Weight(0.0, 0.0)
const Base.abs(w::Weight) = abs(w.d) + abs(w.e) # define abs(Weight)

function partition(order::Int)
    # normal order, G order, W order
    par = [(1, 0, 0),  # order 1
        (2, 0, 0), (1, 1, 0), (1, 0, 1),  #order 2
        (3, 0, 0), (2, 1, 0), (2, 0, 1), (1, 1, 1), (1, 2, 0), (1, 0, 2), #order 3
        (4, 0, 0), (3, 1, 0), (3, 0, 1), (2, 1, 1), (2, 2, 0), (2, 0, 2), (1, 3, 0), (1, 0, 3), (1, 2, 1), (1, 1, 2) #order 4
    ]
    return sort([p for p in par if p[1] + p[2] + p[3] <= order])
end

function neighbor(partitions)
    n = Vector{Tuple{Int,Int}}()
    Nnorm = length(partitions) + 1 # the index of the normalization diagram is the N+1
    for (ip, p) in enumerate(partitions)
        if p[1] == 1 # if there is only one loop, then the diagram can be connected to the normalization diagram
            push!(n, (ip, Nnorm))
        end
        for (idx, np) in enumerate(partitions)
            if idx >= ip
                continue
            end
            if np[1] == p[1] || np[1] == p[1] + 1 || np[1] == p[1] - 1 #the first index is the number of loops
                push!(n, (ip, idx))
            end
        end
    end
    println(n)
    return n
end

include("interaction.jl")

end