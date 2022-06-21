# We work with Rydberg units, length scale Bohr radius a_0, energy scale: Ry
using StaticArrays
using ElectronGas: Parameter
using CompositeGrids
using FeynmanDiagram

include("counterterm.jl")

const D = 2
const beta = 100.0
const rs = 1.0
# const mass2 = 0.3838^2
const mass2 = 1e-6
# const Fs = -0.58545
# const Fs = -0.20633
const Fs = -0.0
const Fa = -0.0
# const massratio = 1.049
const massratio = 1.0
const z = 1.0

const isFock = false

const parafileName = joinpath(@__DIR__, "para.csv") # ROOT/common/para.csv

const interaction = [FeynmanDiagram.Interaction(ChargeCharge, [
    Instant,
    Dynamic
]),]  #instant charge-charge interaction


const paraid = Dict("D" => D, "rs" => rs, "beta" => beta, "mass2" => mass2, "Fs" => Fs,
    "Fa" => Fa, "massratio" => massratio, "z" => z, "isFock" => isFock, "interaction" => interaction
)

const para = Parameter.rydbergUnit(1.0 / beta, rs, D, Λs=mass2)

###### constants ###########
const kF = para.kF
const EF = para.EF
const β = para.β
const me = para.me
const e0 = para.e0
const ϵ0 = para.ϵ0
const dim = para.dim
const μ = para.μ
const NF = para.NF
const spin = para.spin
const maxK = 6 * kF

const qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, maxK], [0.0, 2kF], 16, 0.01 * kF, 8)
const τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)

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
    return sort([p for p in par if p[1] + p[2] + p[3] <= Order])
end


# println("rs=$rs, β=$β, kF=$kF, EF=$EF, mass2=$mass2, NF=$NF, qTF/kF=$(qTF/kF)")
println(para)