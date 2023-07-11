module UEG
using DataStructures
using StaticArrays
using ..Parameters
using ..ElectronGas: Parameter
using ..CompositeGrids

export INL, INR, OUTL, OUTR, DI, EX
export ParaMC, Weight, getK

########### constant that don't needs to be updated frequently ##############
# const Dim = 3
# const Spin = 2
# const IsDynamic = true
# # const IsDynamic = false
# const IsFock = false
############################################################################

# Specify the type of qgrid and τgrid explicitly, otherwise, there will be a type stability issue with interactionDynamic and interactionStatic
# const GridType = CompositeGrids.CompositeG.Composite{Float64,CompositeGrids.SimpleG.Arbitrary{Float64,CompositeGrids.SimpleG.ClosedBound},CompositeGrids.CompositeG.Composite{Float64,CompositeGrids.SimpleG.Log{Float64},CompositeGrids.SimpleG.Uniform{Float64,CompositeGrids.SimpleG.ClosedBound}}}

# NOTE: Reverting to CompositeGrids v0.1.1 type for compatibility with SOSEM JLD2 archives
const GridType = CompositeGrids.CompositeG.Composite{Float64,CompositeGrids.SimpleG.Arbitrary{Float64},CompositeGrids.CompositeG.Composite{Float64,CompositeGrids.SimpleG.Log{Float64},CompositeGrids.SimpleG.Uniform{Float64}}}

@with_kw mutable struct ParaMC
    ### fundamental parameters
    beta::Float64
    rs::Float64
    order::Int = 2
    Fs::Float64 = -0.0
    Fa::Float64 = -0.0
    # δFs = []

    mass2::Float64 = 1e-6
    massratio::Float64 = 1.0

    dim::Int = 3
    spin::Int = 2
    isFock::Bool = false
    isDynamic::Bool = false

    ### MC parameters #######
    # seed::Int = abs(rand(Int)) % 1000000
    # steps::Int = 1e6

    #### derived parameters ###########
    basic::Parameter.Para = Parameter.rydbergUnit(1.0 / beta, rs, dim, Λs=mass2, spin=spin)

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
    qTF::Float64 = basic.qTF

    fs::Float64 = Fs / NFstar
    fa::Float64 = Fa / NFstar

    ##########   effective interaction and counterterm ###############
    qgrid::GridType = CompositeGrid.LogDensedGrid(:uniform, [0.0, maxK], [0.0, 2kF], 16, 0.01 * kF, 8)
    τgrid::GridType = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)

    # ######### only need to be initialized for MC simulation ###########################
    initialized::Bool = false
    dW0::Matrix{Float64} = Matrix{Float64}(undef, length(qgrid), length(τgrid))
    # dW0_f::Matrix{Float64} = Matrix{Float64}(undef, length(qgrid), length(τgrid))
    cRs::Vector{Matrix{Float64}} = []
    # cRs_f::Vector{Matrix{Float64}} = []

    # dW0::Matrix{Float64} = KOdynamic_T(basic, qgrid, τgrid, mass2, massratio, fs, fa)
    # cRs::Vector{Matrix{Float64}} = [counterKO_T(basic, qgrid, τgrid, o, mass2, massratio, fs, fa) for o in 1:order]

    additional = Any[]
end

# struct TestPara
#     e0::Float64
# end

function MCinitialize!(para::ParaMC)
    para.dW0 .= KOdynamic_T(para)
    # para.dW0_f .= KOdynamic_T_df(para)
    for o in 1:para.order-1
        push!(para.cRs, counterKO_T(para; order=o))
        # push!(para.cRs_f, counterKO_T_df(para; order=o))
    end
    para.initialized = true
end

const INL, OUTL, INR, OUTR = 1, 2, 3, 4
const DI, EX = 1, 2

paraid(p::ParaMC) = Dict(
    "order" => p.order,
    "dim" => p.dim,
    "rs" => p.rs,
    "beta" => p.beta,
    "mass2" => p.mass2,
    "Fs" => p.Fs,
    "Fa" => p.Fa,
    "massratio" => p.massratio,
    "spin" => p.spin,
    "isFock" => p.isFock,
    "isDynamic" => p.isDynamic,
)

short(p::ParaMC) = join(["$(k)_$(v)" for (k, v) in sort!(OrderedDict(paraid(p)))], "_")

"""Mapping from ParaMC fields saved in short format to their corresponding types"""
const short_paratypes = Dict(
    "order" => Int,
    "dim" => Int,
    "rs" => Float64,
    "beta" => Float64,
    "mass2" => Float64,
    "Fs" => Float64,
    "Fa" => Float64,
    "massratio" => Float64,
    "spin" => Int,
    "isFock" => Bool,
    "isDynamic" => Bool,
)

"""Constructs a ParaMC from a short string representation."""
function ParaMC(short::String)
    # Split short string into key-value pairs
    short_split = split(short, "_")
    @assert iseven(length(short_split))
    keys = short_split[1:2:end]
    valstrings = short_split[2:2:end]
    # Parse values by type
    vals = []
    for (key, valstring) in zip(keys, valstrings)
        @assert haskey(short_paratypes, key)
        val = parse(short_paratypes[key], valstring)
        push!(vals, val)
    end
    # Construct kwargs (:dim => 3, :rs => 1.0, etc.)
    kwargs = Dict(zip(Symbol.(keys), vals))
    # Construct ParaMC object
    return ParaMC(; kwargs...)
end

function Base.isequal(p1::ParaMC, p2::ParaMC)
    # If either p1 or p2 is initialized, then include initializable
    # fields for purposes of object equality. Otherwise, exclude them.
    if p1.initialized || p2.initialized
        uninitialized_fields = []
    else
        uninitialized_fields = [:dW0, :cRs]
        # uninitialized_fields = [:dW0, dW0_f, :cRs, :cRs_f]
    end
    # p1 and p2 are equal if all fields are equal modulo uninitialized fields
    for field in setdiff(fieldnames(ParaMC), uninitialized_fields)
        if getproperty(p1, field) != getproperty(p1, field)
            return false
        end
        # if field == :qgrid
        #     (getproperty(a, :qgrid) ≈ getproperty(b, :qgrid)) == false && return false
        # else
        #     getproperty(a, field) != getproperty(b, field) && return false
        # end
    end
    return true
end
Base.:(==)(p1::ParaMC, p2::ParaMC) = Base.isequal(p1, p2)

"""
Hard-coded counterterm partitions for diagrams of max order `order` and minimal
loop order `offset`, given in the form (n_loop, n_μ, n_λ). The default offset
corresponds to partitions of the self-energy, where the minimal loop order is 1.
"""
function partition(order::Int; offset::Int=1)
    # normal order, G order, W order
    par = [
        # order 0
        (0, 0, 0),
        # order 1
        (1, 0, 0), (0, 1, 0), (0, 0, 1),
        # order 2
        (2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1), (0, 2, 0), (0, 0, 2),
        # order 3
        (3, 0, 0), (2, 1, 0), (2, 0, 1), (1, 1, 1), (1, 2, 0),
        (1, 0, 2), (0, 3, 0), (0, 0, 3), (0, 2, 1), (0, 1, 2),
        #order 4
        (4, 0, 0), (3, 1, 0), (3, 0, 1), (2, 2, 0), (2, 1, 1), (2, 0, 2), (1, 3, 0), (1, 2, 1),
        (1, 1, 2), (1, 0, 3), (0, 4, 0), (0, 3, 1), (0, 2, 2), (0, 1, 3), (0, 0, 4),
        #order 5
        (5, 0, 0), (4, 1, 0), (4, 0, 1), (3, 2, 0), (3, 1, 1), (3, 0, 2), (2, 3, 0), (2, 2, 1),
        (2, 1, 2), (2, 0, 3), (1, 4, 0), (1, 3, 1), (1, 2, 2), (1, 1, 3), (1, 0, 4), (0, 5, 0),
        (0, 4, 1), (0, 3, 2), (0, 2, 3), (0, 1, 4), (0, 0, 5),
    ]
    return sort([(p[1] + offset, p[2], p[3]) for p in par if p[1] + p[2] + p[3] <= order - offset])
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

function getK(amp, dim::Int, idx::Int=1)
    @assert 1 <= idx <= dim
    k = zeros(dim)
    k[idx] = amp
    return k
end

include("interaction.jl")

end