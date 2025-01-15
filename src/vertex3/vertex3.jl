module Ver3

using JLD2, CSV
using Printf, LinearAlgebra, DataFrames
using ..StaticArrays
using ..Parameters
using ..CompositeGrids
using ..ElectronGas
using ..MCIntegration
using ..Lehmann

using ..FeynmanDiagram
import ..FeynmanDiagram.FrontEnds: TwoBodyChannel, Alli, PHr, PHEr, PPr, AnyChan
import ..FeynmanDiagram.FrontEnds: Filter, NoHartree, NoFock, DirectOnly, Wirreducible, Girreducible, NoBubble, Proper
import ..FeynmanDiagram.FrontEnds: Response, Composite, ChargeCharge, SpinSpin, UpUp, UpDown
import ..FeynmanDiagram.FrontEnds: AnalyticProperty, Instant, Dynamic
import ..FeynmanDiagram.Parquet: DiagPara, Ver4Diag, PolarDiag, Ver3Diag
using ..Measurements
using ..Diagram

using ..UEG
using ..Propagator
import ..Propagator: LeafStateADDynamic

# import ..ExprTreeF64
import ..Weight

function diagPara(para::ParaMC, order, filter, transferLoop)
    inter = [FeynmanDiagram.Interaction(ChargeCharge, para.isDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    return DiagParaF64(
        type=Ver3Diag,
        innerLoopNum=order,
        hasTau=true,
        spin=para.spin,
        firstLoopIdx=3,
        interaction=inter,
        filter=filter,
        transferLoop=transferLoop
    )
end

function diagram(paramc::ParaMC, _partition::Vector{T};
    filter=[
        NoHartree,
        # Girreducible,
        Proper,   #one interaction irreduble diagrams or not
        # NoBubble, #allow the bubble diagram or not
    ]
) where {T}
    # println("Build the vertex4 diagrams into an expression tree ...")
    # _partition = UEG.partition(order)
    # println("Diagram set: ", _partition)

    dim = paramc.dim
    Kin, Qout = zeros(16), zeros(16)
    Qout[1], Kin[2] = 1.0, 1.0
    legK = [Qout, Kin]

    diag = Vector{ExprTreeF64}()
    diagpara = Vector{DiagParaF64}()
    partition = Vector{T}()
    for p in _partition
        para = diagPara(paramc, p[1], filter, Qout)
        d::Vector{Diagram{Float64}} = Parquet.vertex3(para, legK).diagram
        d = DiagTree.derivative(d, BareGreenId, p[2], index=1)
        d = DiagTree.derivative(d, BareInteractionId, p[3], index=2)

        # the Taylor expansion should be d^n f(x) / dx^n / n!, so there is a factor of 1/n! for each derivative
        for _d in d
            _d.factor *= 1 / factorial(p[2]) / factorial(p[3])
        end
        if isempty(d) == false
            if paramc.isFock # remove the Fock subdiagrams
                DiagTree.removeHartreeFock!(d)
            end
            push!(diagpara, para)
            push!(diag, ExprTree.build(d, dim))
            push!(partition, p)
        else
            @warn("partition $p doesn't have any diagram. It will be ignored.")
        end
    end

    # diag = [ExprTree.build(d) for d in ver4]    #expression tree representation of diagrams 
    rootuu = [[idx for idx in d.root if d.node.object[idx].para.response == UpUp] for d in diag] #select the diagram with upup
    rootud = [[idx for idx in d.root if d.node.object[idx].para.response == UpDown] for d in diag] #select the diagram with updown
    #assign the external Tau to the corresponding diagrams
    extTuu = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootuu)]
    extTud = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootud)]
    return (partition, diagpara, diag, [rootuu, rootud], [extTuu, extTud])
end

@inline function phaseF(varT, extT, nin, nout, β)
    # println(extT)
    tb, tfin, tfout = varT[extT[1]], varT[extT[2]], varT[extT[3]]
    win, wout = π * (2nin + 1) / β, π * (2nout + 1) / β
    wq = win - wout
    return cos(tfin * win - tfout * wout - tb * wq)

    # if (idx == 1)
    #     return cos(π / β * ((2tb) - (tfin + tfout)))
    #     # return cos(π / β * ((2tb) - 3 * tfin + tfout))
    #     # return cos(π / β * ((2tb) - 3 * tfin + tfout))
    # else
    #     return cos(π / β * (tfin - tfout))
    # end
end

@inline function phaseC(varT, extT, nin, nout, β)
    # println(extT)
    tb, tfin, tfout = varT[extT[1]], varT[extT[2]], varT[extT[3]]
    win, wout = π * (2nin + 1) / β, π * (2nout + 1) / β
    wq = win - wout
    return exp(-1im * (tfin * win - tfout * wout - tb * wq))
end

function diagram_loadinfo(paramc::ParaMC, _partition::Vector{T};
    filter=[NoHartree], transferLoop=nothing,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"), filename="extvars_vertex3.jld2"
) where {T}
    diagpara = Vector{DiagPara}()
    extT_labels = Vector{Vector{Int}}[]
    spin_conventions = Vector{Response}[]

    fname = joinpath(root_dir, filename)
    jldopen(fname, "r") do f
        for p in _partition
            key_str = join(string.(p))
            if key_str in keys(f)
                extT, ext_spin = f[key_str]
                push!(diagpara, Diagram.diagPara(Ver3Diag, paramc.isDynamic, p[1], paramc.spin, filter, transferLoop))
                push!(extT_labels, extT)
                push!(spin_conventions, ext_spin)
            else
                error("$(key_str) not found in $(fname)")
            end
        end
    end
    return (_partition, diagpara, extT_labels, spin_conventions)
end
# @inline function phase(varT, extT, n, β)
#     # println(extT)
#     return phase(varT, extT, n[1], n[2], n[3], β)
# end

# include("ver3KW.jl")
include("ver3KW_Clib.jl")
# include("ver3angleavg.jl")
include("source_codeParquetAD/Cwrapper_vertex3_ParquetAD.jl")


const evalfuncParquetAD_vertex3_map = Dict(
    (1, 0, 0) => eval_vertex3_ParquetAD100!,
    (1, 0, 1) => eval_vertex3_ParquetAD101!,
    (1, 0, 2) => eval_vertex3_ParquetAD102!,
    (1, 0, 3) => eval_vertex3_ParquetAD103!,
    (1, 0, 4) => eval_vertex3_ParquetAD104!,
    (1, 0, 5) => eval_vertex3_ParquetAD105!,
    (1, 1, 0) => eval_vertex3_ParquetAD110!,
    (1, 1, 1) => eval_vertex3_ParquetAD111!,
    (1, 1, 2) => eval_vertex3_ParquetAD112!,
    (1, 1, 3) => eval_vertex3_ParquetAD113!,
    (1, 1, 4) => eval_vertex3_ParquetAD114!,
    (1, 2, 0) => eval_vertex3_ParquetAD120!,
    (1, 2, 1) => eval_vertex3_ParquetAD121!,
    (1, 2, 2) => eval_vertex3_ParquetAD122!,
    (1, 2, 3) => eval_vertex3_ParquetAD123!,
    (1, 3, 0) => eval_vertex3_ParquetAD130!,
    (1, 3, 1) => eval_vertex3_ParquetAD131!,
    (1, 3, 2) => eval_vertex3_ParquetAD132!,
    (1, 4, 0) => eval_vertex3_ParquetAD140!,
    (1, 4, 1) => eval_vertex3_ParquetAD141!,
    (1, 5, 0) => eval_vertex3_ParquetAD150!,
    (2, 0, 0) => eval_vertex3_ParquetAD200!,
    (2, 0, 1) => eval_vertex3_ParquetAD201!,
    (2, 0, 2) => eval_vertex3_ParquetAD202!,
    (2, 0, 3) => eval_vertex3_ParquetAD203!,
    (2, 0, 4) => eval_vertex3_ParquetAD204!,
    (2, 1, 0) => eval_vertex3_ParquetAD210!,
    (2, 1, 1) => eval_vertex3_ParquetAD211!,
    (2, 1, 2) => eval_vertex3_ParquetAD212!,
    (2, 1, 3) => eval_vertex3_ParquetAD213!,
    (2, 2, 0) => eval_vertex3_ParquetAD220!,
    (2, 2, 1) => eval_vertex3_ParquetAD221!,
    (2, 2, 2) => eval_vertex3_ParquetAD222!,
    (2, 3, 0) => eval_vertex3_ParquetAD230!,
    (2, 3, 1) => eval_vertex3_ParquetAD231!,
    (2, 4, 0) => eval_vertex3_ParquetAD240!,
    (3, 0, 0) => eval_vertex3_ParquetAD300!,
    (3, 0, 1) => eval_vertex3_ParquetAD301!,
    (3, 0, 2) => eval_vertex3_ParquetAD302!,
    (3, 0, 3) => eval_vertex3_ParquetAD303!,
    (3, 1, 0) => eval_vertex3_ParquetAD310!,
    (3, 1, 1) => eval_vertex3_ParquetAD311!,
    (3, 1, 2) => eval_vertex3_ParquetAD312!,
    (3, 2, 0) => eval_vertex3_ParquetAD320!,
    (3, 2, 1) => eval_vertex3_ParquetAD321!,
    (3, 3, 0) => eval_vertex3_ParquetAD330!,
    (4, 0, 0) => eval_vertex3_ParquetAD400!,
    (4, 0, 1) => eval_vertex3_ParquetAD401!,
    (4, 0, 2) => eval_vertex3_ParquetAD402!,
    (4, 1, 0) => eval_vertex3_ParquetAD410!,
    (4, 1, 1) => eval_vertex3_ParquetAD411!,
    (4, 2, 0) => eval_vertex3_ParquetAD420!,
    (5, 0, 0) => eval_vertex3_ParquetAD500!,
    (5, 0, 1) => eval_vertex3_ParquetAD501!,
    (5, 1, 0) => eval_vertex3_ParquetAD510!,
    (6, 0, 0) => eval_vertex3_ParquetAD600!
)

end