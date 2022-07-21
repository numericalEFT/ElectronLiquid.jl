module UEG
using StaticArrays
using Parameters
using ElectronGas: Parameter
using CompositeGrids

export INL, INR, OUTL, OUTR, DI, EX
export ParaMC, Weight

########### constant that don't needs to be updated frequently ##############
const Dim = 3
const Spin = 2
const IsDynamic = true
const IsFock = false
############################################################################

#specify the type of qgrid and τgrid explicitly, otherwise, there will be a type stability issue with interactionDynamic and interactionStatic
const GridType = CompositeGrids.CompositeG.Composite{Float64,CompositeGrids.SimpleG.Arbitrary{Float64},CompositeGrids.CompositeG.Composite{Float64,CompositeGrids.SimpleG.Log{Float64},CompositeGrids.SimpleG.Uniform{Float64}}}

@with_kw struct ParaMC
    ### fundamental parameters
    beta::Float64
    rs::Float64
    order::Int = 2
    Fs::Float64 = -0.0
    Fa::Float64 = -0.0
    # δFs = []

    mass2::Float64 = 1e-6
    massratio::Float64 = 1.0

    ### MC parameters #######
    # seed::Int = abs(rand(Int)) % 1000000
    # steps::Int = 1e6

    #### derived parameters ###########
    basic::Parameter.Para = Parameter.rydbergUnit(1.0 / beta, rs, Dim, Λs=mass2, spin=Spin)

    dim::Int = basic.dim
    spin::Int = basic.spin
    kF::Float64 = basic.kF
    EF::Float64 = basic.EF
    β::Float64 = basic.β
    maxK::Float64 = 6 * basic.kF
    me::Float64 = basic.me
    ϵ0::Float64 = basic.ϵ0
    e0::Float64 = basic.e0
    μ::Float64 = basic.μ
    NF::Float64 = basic.NF
    NFstar::Float64 = basic.NF * massratio

    fs::Float64 = Fs / NFstar
    fa::Float64 = Fa / NFstar

    ##########   effective interaction and counterterm ###############
    qgrid::GridType = CompositeGrid.LogDensedGrid(:uniform, [0.0, maxK], [0.0, 2kF], 16, 0.01 * kF, 8)
    τgrid::GridType = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)

    ######### only need to be initialized for MC simulation ###########################
    dW0::Matrix{Float64} = Matrix{Float64}(undef, length(qgrid), length(τgrid))
    cRs::Vector{Matrix{Float64}} = []

    # dW0::Matrix{Float64} = KOdynamic_T(basic, qgrid, τgrid, mass2, massratio, fs, fa)
    # cRs::Vector{Matrix{Float64}} = [counterKO_T(basic, qgrid, τgrid, o, mass2, massratio, fs, fa) for o in 1:order]

    additional = Any[]
end

function MCinitialize!(para::ParaMC)
    para.dW0 .= KOdynamic_T(para)
    for o in 1:para.order-1
        push!(para.cRs, counterKO_T(para; order=o))
    end
end

const INL, OUTL, INR, OUTR = 1, 2, 3, 4
const DI, EX = 1, 2

paraid(p::ParaMC) = Dict(
    "dim" => p.dim,
    "rs" => p.rs,
    "beta" => p.beta,
    "mass2" => p.mass2,
    "Fs" => p.Fs,
    "Fa" => p.Fa,
    "massratio" => p.massratio,
    "spin" => p.spin,
    "isFock" => IsFock,
    "isDynamic" => IsDynamic,
)

short(p::ParaMC) = join(["$(k)_$(v)" for (k, v) in sort(paraid(p))], "_")


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