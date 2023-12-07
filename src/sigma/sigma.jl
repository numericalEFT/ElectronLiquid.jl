module Sigma
using Cuba
using JLD2, CSV

using Printf, LinearAlgebra, DataFrames
using ..CompositeGrids
using ..ElectronGas
using ..MCIntegration
using ..Lehmann

using ..FeynmanDiagram
using ..Measurements

using ..UEG
using ..Propagator
import ..ExprTreeF64

function diagPara(para::ParaMC, order::Int, filter)
    inter = [FeynmanDiagram.Interaction(ChargeCharge, para.isDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    DiagParaF64(
        type=SigmaDiag,
        innerLoopNum=order,
        hasTau=true,
        # loopDim=para.dim,
        spin=para.spin,
        firstLoopIdx=2,
        interaction=inter,
        filter=filter
    )
end

function diagram(paramc::ParaMC, _partition::Vector{T};
    filter=[
        FeynmanDiagram.NoHartree,
        # Girreducible,
        # Proper,   #one interaction irreduble diagrams or not
        # NoBubble, #allow the bubble diagram or not
    ],
    dR=false # whether to use the derivative of the renormalized interaction, for the RG purpose
) where {T}
    println("Build the sigma diagrams into an expression tree ...")
    println("Diagram set: ", _partition)

    dim = paramc.dim
    diag = Vector{ExprTreeF64}()
    diagpara = Vector{DiagParaF64}()
    partition = Vector{T}()
    diagrams = Vector{Diagram{Float64}}()
    for p in _partition
        para = diagPara(paramc, p[1], filter)
        sd::Vector{Diagram{Float64}} = Parquet.sigma(para).diagram
        sdp = DiagTree.derivative(sd, BareGreenId, p[2], index=1)
        sdpp = DiagTree.derivative(sdp, BareInteractionId, p[3], index=2)
        if dR
            sdpp = DiagTree.derivative(sdpp, BareInteractionId, 1, index=3)
        end
        # the Taylor expansion should be d^n f(x) / dx^n / n!, so there is a factor of 1/n! for each derivative
        for d in sdpp
            d.factor *= 1 / factorial(p[2]) / factorial(p[3])
        end
        if isempty(sdpp) == false
            if paramc.isFock && (p != (1, 0, 0)) # the Fock diagram itself should not be removed
                DiagTree.removeHartreeFock!(sdpp)
            end
            push!(diagpara, para)
            push!(partition, p)
            push!(diag, ExprTree.build(sdpp, dim))
            append!(diagrams, sdpp)
        else
            @warn("partition $p doesn't have any diagram. It will be ignored.")
        end
    end

    # @time diag = [ExprTree.build(diags) for diags in sigma] # DiagTree to ExprTree
    root = [d.root for d in diag] #get the list of root nodes
    #assign the external Tau to the corresponding diagrams
    extT = [[diag[ri].node.object[idx].para.extT::Tuple{Int,Int} for idx in r] for (ri, r) in enumerate(root)]
    # result = (partition, diagpara, diag, root, extT)
    result = (partition, diagpara, diag, root, extT, diagrams)
    return result
end

function diagramParquet(paramc::ParaMC, _partition::Vector{T}; filter=[FeynmanDiagram.NoHartree]) where {T}
    diagpara = Vector{DiagParaF64}()
    extT_labels = Vector{Vector{Int}}[]

    jldopen(joinpath(@__DIR__, "source_codeGV", "extT_ParquetAD.jld2"), "r") do f
        for p in _partition
            push!(diagpara, diagPara(paramc, p[1], filter))
            key_str = join(string.(p))
            push!(extT_labels, f[key_str])
        end
    end
    return (_partition, diagpara, extT_labels)
end

function diagramParquet(paramc::ParaMC, _partition::Vector{T}, spinPolarPara::Float64;
    filter=[FeynmanDiagram.NoHartree]) where {T}
    diagpara = Vector{DiagParaF64}()
    FeynGraphs = FeynmanDiagram.diagdict_parquet(:sigma, _partition, spinPolarPara=spinPolarPara, filter=filter, isDynamic=paramc.isDynamic)
    extT_labels = Vector{Vector{Int}}[]
    for p in _partition
        push!(diagpara, diagPara(paramc, p[1], filter))
        push!(extT_labels, FeynGraphs[p][2])
    end
    return (_partition, diagpara, FeynGraphs, extT_labels)
end

function diagramGV(paramc::ParaMC, _partition::Vector{T}; filter=[FeynmanDiagram.NoHartree]) where {T}
    diagpara = Vector{DiagParaF64}()
    extT_labels = Vector{Vector{Int}}[]

    for p in _partition
        push!(diagpara, diagPara(paramc, p[1], filter))
        if p[1] == 1
            push!(extT_labels, [[1, 1]])
        else
            push!(extT_labels, [[1, 1], [1, 2]])
        end
    end

    return (_partition, diagpara, extT_labels)
end

function diagramGV(paramc::ParaMC, _partition::Vector{T}, spinPolarPara::Float64;
    filter=[FeynmanDiagram.NoHartree]) where {T}
    diagpara = Vector{DiagParaF64}()
    extT_labels = Vector{Vector{Int}}[]
    for p in _partition
        push!(diagpara, diagPara(paramc, p[1], filter))
        if p[1] == 1
            push!(extT_labels, [[1, 1]])
        else
            push!(extT_labels, [[1, 1], [1, 2]])
        end
    end
    FeynGraphs, labelProd = FeynmanDiagram.diagdictGV(:sigma, _partition, spinPolarPara=spinPolarPara)

    return (_partition, diagpara, FeynGraphs, labelProd, extT_labels)
end

@inline function phase(varT, extT, l, β)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return exp(1im * π * (2l + 1) / β * (tout - tin))
end

#include("sigma_generic.jl")
include("sigma_dk.jl")
include("sigmaKW.jl")
include("sigmaCuba.jl")
include("sigmaVegas.jl")
include("sigmaGV.jl")
include("sigmaGV_AD.jl")
include("sigmaGV_compile.jl")


function MC(para; kgrid=[para.kF,], ngrid=[-1, 0, 1], neval=1e6, reweight_goal=nothing,
    spinPolarPara::Float64=0.0, # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
    filename::Union{String,Nothing}=nothing, partition=UEG.partition(para.order), diagtype=:Parquet,
    isLayered2D=false, # whether to use the screened Coulomb interaction in 2D or not 
    isClib=true # whether to use compiled C library to calculate the Feynman diagram weight or not (only for spin-unpolarized case now)
)
    kF = para.kF
    neighbor = UEG.neighbor(partition)

    if isLayered2D
        @assert (para.dim == 2) && diagtype == :GV "Only 2D and GV diagrams supports the tanh screened Coulomb interaction"
    end

    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            reweight_factor = 2.0^(2order + sOrder + vOrder - 2)
            if (order, sOrder, vOrder) == (1, 0, 0)
                reweight_factor = 4.0
            end
            push!(reweight_goal, reweight_factor)
        end
        push!(reweight_goal, 4.0)
    end

    if diagtype == :GV
        if spinPolarPara == 0.0 && isClib
            diagram = Sigma.diagramGV(para, partition)
            sigma, result = Sigma.GV_Clib(para, diagram;
                isLayered2D=isLayered2D,
                neighbor=neighbor, reweight_goal=reweight_goal,
                kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread)
        else
            diagram = Sigma.diagramGV(para, partition, spinPolarPara)
            sigma, result = Sigma.GV(para, diagram;
                isLayered2D=isLayered2D,
                neighbor=neighbor, reweight_goal=reweight_goal,
                kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread)
        end
    elseif diagtype == :Parquet
        if spinPolarPara == 0.0 && isClib
            diagram = Sigma.diagramParquet(para, partition)
            sigma, result = Sigma.ParquetAD_Clib(para, diagram;
                isLayered2D=isLayered2D,
                neighbor=neighbor, reweight_goal=reweight_goal,
                kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread)
        else
            # diagram = Sigma.diagram(para, partition)
            # sigma, result = Sigma.KW(para, diagram;
            #     neighbor=neighbor, reweight_goal=reweight_goal,
            #     kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread)
            diagram = Sigma.diagramParquet(para, partition, spinPolarPara)
            sigma, result = Sigma.ParquetAD(para, diagram;
                isLayered2D=isLayered2D,
                neighbor=neighbor, reweight_goal=reweight_goal,
                kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread)
        end
    else
        error("unknown diagrams' generated type")
    end

    if isnothing(sigma) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (ngrid, kgrid, sigma)
            end
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            @printf("%10s  %10s   %10s   %10s   %10s \n", "q/kF", "real(avg)", "err", "imag(avg)", "err")
            r, i = real(sigma[key]), imag(sigma[key])
            for (in, n) in enumerate(ngrid)
                println("n = $n")
                for (iq, q) in enumerate(kgrid)
                    @printf("%10.6f  %10.6f ± %10.6f   %10.6f ± %10.6f\n", q[1] / kF, r[in, iq].val, r[in, iq].err, i[in, iq].val, i[in, iq].err)
                end
            end
        end
    end
    return sigma, result
end

include("source_codeGV/Cwrapper_sigmaGV.jl")
include("source_codeGV/Cwrapper_sigmaParquetAD.jl")
include("source_codeGV/Cwrapper_sigmaParquetAD_dk.jl")
const evalfuncGV_map = Dict(
    (1, 0, 0) => eval_sigmaGV100!,
    (1, 0, 1) => eval_sigmaGV101!,
    (1, 0, 2) => eval_sigmaGV102!,
    (1, 0, 3) => eval_sigmaGV103!,
    (1, 0, 4) => eval_sigmaGV104!,
    (1, 0, 5) => eval_sigmaGV105!,
    (1, 1, 0) => eval_sigmaGV110!,
    (1, 1, 1) => eval_sigmaGV111!,
    (1, 1, 2) => eval_sigmaGV112!,
    (1, 1, 3) => eval_sigmaGV113!,
    (1, 1, 4) => eval_sigmaGV114!,
    (1, 2, 0) => eval_sigmaGV120!,
    (1, 2, 1) => eval_sigmaGV121!,
    (1, 2, 2) => eval_sigmaGV122!,
    (1, 2, 3) => eval_sigmaGV123!,
    (1, 3, 0) => eval_sigmaGV130!,
    (1, 3, 1) => eval_sigmaGV131!,
    (1, 3, 2) => eval_sigmaGV132!,
    (1, 4, 0) => eval_sigmaGV140!,
    (1, 4, 1) => eval_sigmaGV141!,
    (1, 5, 0) => eval_sigmaGV150!,
    (2, 0, 0) => eval_sigmaGV200!,
    (2, 0, 1) => eval_sigmaGV201!,
    (2, 0, 2) => eval_sigmaGV202!,
    (2, 0, 3) => eval_sigmaGV203!,
    (2, 0, 4) => eval_sigmaGV204!,
    (2, 1, 0) => eval_sigmaGV210!,
    (2, 1, 1) => eval_sigmaGV211!,
    (2, 1, 2) => eval_sigmaGV212!,
    (2, 1, 3) => eval_sigmaGV213!,
    (2, 2, 0) => eval_sigmaGV220!,
    (2, 2, 1) => eval_sigmaGV221!,
    (2, 2, 2) => eval_sigmaGV222!,
    (2, 3, 0) => eval_sigmaGV230!,
    (2, 3, 1) => eval_sigmaGV231!,
    (2, 4, 0) => eval_sigmaGV240!,
    (3, 0, 0) => eval_sigmaGV300!,
    (3, 0, 1) => eval_sigmaGV301!,
    (3, 0, 2) => eval_sigmaGV302!,
    (3, 0, 3) => eval_sigmaGV303!,
    (3, 1, 0) => eval_sigmaGV310!,
    (3, 1, 1) => eval_sigmaGV311!,
    (3, 1, 2) => eval_sigmaGV312!,
    (3, 2, 0) => eval_sigmaGV320!,
    (3, 2, 1) => eval_sigmaGV321!,
    (3, 3, 0) => eval_sigmaGV330!,
    (4, 0, 0) => eval_sigmaGV400!,
    (4, 0, 1) => eval_sigmaGV401!,
    (4, 0, 2) => eval_sigmaGV402!,
    (4, 1, 0) => eval_sigmaGV410!,
    (4, 1, 1) => eval_sigmaGV411!,
    (4, 2, 0) => eval_sigmaGV420!,
    (5, 0, 0) => eval_sigmaGV500!,
    (5, 0, 1) => eval_sigmaGV501!,
    (5, 1, 0) => eval_sigmaGV510!,
    (6, 0, 0) => eval_sigmaGV600!
)
const evalfuncParquetAD_map = Dict(
    (1, 0, 0) => eval_sigmaParquetAD100!,
    (1, 0, 1) => eval_sigmaParquetAD101!,
    (1, 0, 2) => eval_sigmaParquetAD102!,
    (1, 0, 3) => eval_sigmaParquetAD103!,
    (1, 0, 4) => eval_sigmaParquetAD104!,
    (1, 0, 5) => eval_sigmaParquetAD105!,
    (1, 1, 0) => eval_sigmaParquetAD110!,
    (1, 1, 1) => eval_sigmaParquetAD111!,
    (1, 1, 2) => eval_sigmaParquetAD112!,
    (1, 1, 3) => eval_sigmaParquetAD113!,
    (1, 1, 4) => eval_sigmaParquetAD114!,
    (1, 2, 0) => eval_sigmaParquetAD120!,
    (1, 2, 1) => eval_sigmaParquetAD121!,
    (1, 2, 2) => eval_sigmaParquetAD122!,
    (1, 2, 3) => eval_sigmaParquetAD123!,
    (1, 3, 0) => eval_sigmaParquetAD130!,
    (1, 3, 1) => eval_sigmaParquetAD131!,
    (1, 3, 2) => eval_sigmaParquetAD132!,
    (1, 4, 0) => eval_sigmaParquetAD140!,
    (1, 4, 1) => eval_sigmaParquetAD141!,
    (1, 5, 0) => eval_sigmaParquetAD150!,
    (2, 0, 0) => eval_sigmaParquetAD200!,
    (2, 0, 1) => eval_sigmaParquetAD201!,
    (2, 0, 2) => eval_sigmaParquetAD202!,
    (2, 0, 3) => eval_sigmaParquetAD203!,
    (2, 0, 4) => eval_sigmaParquetAD204!,
    (2, 1, 0) => eval_sigmaParquetAD210!,
    (2, 1, 1) => eval_sigmaParquetAD211!,
    (2, 1, 2) => eval_sigmaParquetAD212!,
    (2, 1, 3) => eval_sigmaParquetAD213!,
    (2, 2, 0) => eval_sigmaParquetAD220!,
    (2, 2, 1) => eval_sigmaParquetAD221!,
    (2, 2, 2) => eval_sigmaParquetAD222!,
    (2, 3, 0) => eval_sigmaParquetAD230!,
    (2, 3, 1) => eval_sigmaParquetAD231!,
    (2, 4, 0) => eval_sigmaParquetAD240!,
    (3, 0, 0) => eval_sigmaParquetAD300!,
    (3, 0, 1) => eval_sigmaParquetAD301!,
    (3, 0, 2) => eval_sigmaParquetAD302!,
    (3, 0, 3) => eval_sigmaParquetAD303!,
    (3, 1, 0) => eval_sigmaParquetAD310!,
    (3, 1, 1) => eval_sigmaParquetAD311!,
    (3, 1, 2) => eval_sigmaParquetAD312!,
    (3, 2, 0) => eval_sigmaParquetAD320!,
    (3, 2, 1) => eval_sigmaParquetAD321!,
    (3, 3, 0) => eval_sigmaParquetAD330!,
    (4, 0, 0) => eval_sigmaParquetAD400!,
    (4, 0, 1) => eval_sigmaParquetAD401!,
    (4, 0, 2) => eval_sigmaParquetAD402!,
    (4, 1, 0) => eval_sigmaParquetAD410!,
    (4, 1, 1) => eval_sigmaParquetAD411!,
    (4, 2, 0) => eval_sigmaParquetAD420!,
    (5, 0, 0) => eval_sigmaParquetAD500!,
    (5, 0, 1) => eval_sigmaParquetAD501!,
    (5, 1, 0) => eval_sigmaParquetAD510!,
    (6, 0, 0) => eval_sigmaParquetAD600!
)
const evalfuncParquetAD_dk_map = Dict(
    (1, 0, 0) => eval_sigmaParquetAD_dk100!,
    (1, 0, 1) => eval_sigmaParquetAD_dk101!,
    (1, 0, 2) => eval_sigmaParquetAD_dk102!,
    (1, 0, 3) => eval_sigmaParquetAD_dk103!,
    (1, 0, 4) => eval_sigmaParquetAD_dk104!,
    (1, 0, 5) => eval_sigmaParquetAD_dk105!,
    (1, 1, 0) => eval_sigmaParquetAD_dk110!,
    (1, 1, 1) => eval_sigmaParquetAD_dk111!,
    (1, 1, 2) => eval_sigmaParquetAD_dk112!,
    (1, 1, 3) => eval_sigmaParquetAD_dk113!,
    (1, 1, 4) => eval_sigmaParquetAD_dk114!,
    (1, 2, 0) => eval_sigmaParquetAD_dk120!,
    (1, 2, 1) => eval_sigmaParquetAD_dk121!,
    (1, 2, 2) => eval_sigmaParquetAD_dk122!,
    (1, 2, 3) => eval_sigmaParquetAD_dk123!,
    (1, 3, 0) => eval_sigmaParquetAD_dk130!,
    (1, 3, 1) => eval_sigmaParquetAD_dk131!,
    (1, 3, 2) => eval_sigmaParquetAD_dk132!,
    (1, 4, 0) => eval_sigmaParquetAD_dk140!,
    (1, 4, 1) => eval_sigmaParquetAD_dk141!,
    (1, 5, 0) => eval_sigmaParquetAD_dk150!,
    (2, 0, 0) => eval_sigmaParquetAD_dk200!,
    (2, 0, 1) => eval_sigmaParquetAD_dk201!,
    (2, 0, 2) => eval_sigmaParquetAD_dk202!,
    (2, 0, 3) => eval_sigmaParquetAD_dk203!,
    (2, 0, 4) => eval_sigmaParquetAD_dk204!,
    (2, 1, 0) => eval_sigmaParquetAD_dk210!,
    (2, 1, 1) => eval_sigmaParquetAD_dk211!,
    (2, 1, 2) => eval_sigmaParquetAD_dk212!,
    (2, 1, 3) => eval_sigmaParquetAD_dk213!,
    (2, 2, 0) => eval_sigmaParquetAD_dk220!,
    (2, 2, 1) => eval_sigmaParquetAD_dk221!,
    (2, 2, 2) => eval_sigmaParquetAD_dk222!,
    (2, 3, 0) => eval_sigmaParquetAD_dk230!,
    (2, 3, 1) => eval_sigmaParquetAD_dk231!,
    (2, 4, 0) => eval_sigmaParquetAD_dk240!,
    (3, 0, 0) => eval_sigmaParquetAD_dk300!,
    (3, 0, 1) => eval_sigmaParquetAD_dk301!,
    (3, 0, 2) => eval_sigmaParquetAD_dk302!,
    (3, 0, 3) => eval_sigmaParquetAD_dk303!,
    (3, 1, 0) => eval_sigmaParquetAD_dk310!,
    (3, 1, 1) => eval_sigmaParquetAD_dk311!,
    (3, 1, 2) => eval_sigmaParquetAD_dk312!,
    (3, 2, 0) => eval_sigmaParquetAD_dk320!,
    (3, 2, 1) => eval_sigmaParquetAD_dk321!,
    (3, 3, 0) => eval_sigmaParquetAD_dk330!,
    (4, 0, 0) => eval_sigmaParquetAD_dk400!,
    (4, 0, 1) => eval_sigmaParquetAD_dk401!,
    (4, 0, 2) => eval_sigmaParquetAD_dk402!,
    (4, 1, 0) => eval_sigmaParquetAD_dk410!,
    (4, 1, 1) => eval_sigmaParquetAD_dk411!,
    (4, 2, 0) => eval_sigmaParquetAD_dk420!,
    (5, 0, 0) => eval_sigmaParquetAD_dk500!,
    (5, 0, 1) => eval_sigmaParquetAD_dk501!,
    (5, 1, 0) => eval_sigmaParquetAD_dk510!,
    (6, 0, 0) => eval_sigmaParquetAD_dk600!
)
end
