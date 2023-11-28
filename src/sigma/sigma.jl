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
        loopDim=para.dim,
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

    diag = Vector{ExprTreeF64}()
    diagpara = Vector{DiagParaF64}()
    partition = Vector{T}()
    # diagrams = Vector{Diagram{Float64}}()
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
            push!(diag, ExprTree.build(sdpp))
            # append!(diagrams, sdpp)
        else
            @warn("partition $p doesn't have any diagram. It will be ignored.")
        end
    end

    # @time diag = [ExprTree.build(diags) for diags in sigma] # DiagTree to ExprTree
    root = [d.root for d in diag] #get the list of root nodes
    #assign the external Tau to the corresponding diagrams
    extT = [[diag[ri].node.object[idx].para.extT::Tuple{Int,Int} for idx in r] for (ri, r) in enumerate(root)]
    result = (partition, diagpara, diag, root, extT)
    # result = (partition, diagpara, diag, root, extT, diagrams)
    return result
end

function diagramGV(paramc::ParaMC, _partition::Vector{T}; 
    filter=[FeynmanDiagram.NoHartree], root_dir=joinpath(@__DIR__, "source_codeGV/")) where {T}
    diagpara = Vector{DiagParaF64}()
    extT_labels = Vector{Vector{Int}}[]
    for p in _partition
        push!(diagpara, diagPara(paramc, p[1], filter))
        if p[1] == 1
            push!(extT_labels, [[1,1]])
        else
            push!(extT_labels, [[1,1], [1,2]])
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
            push!(extT_labels, [[1,1]])
        else
            push!(extT_labels, [[1,1], [1,2]])
        end
    end
    FeynGraphs, labelProd = FeynmanDiagram.diagdictGV(:sigma, _partition, spinPolarPara=spinPolarPara)
    
    return (_partition, diagpara, FeynGraphs, labelProd, extT_labels)
end

function GVcompile_toFile(maxOrder::Int, FeynGraphs, labelProd, root_dir=joinpath(@__DIR__, "source_codeGV/"))
    leaf_maps = Vector{Dict{Int,FeynmanGraph}}()
    ### compile and save the generated Julia function to a source code file.
    for order in 1:maxOrder
        partition = UEG.partition_order(order)
        for key in partition
            key_str = join(string.(key))
            leafmap = Compilers.compile_code(FeynGraphs[key][1], root_dir * "func_sigmaGV_o$order.jl"; func_name="eval_graph$(key_str)!")
            push!(leaf_maps, leafmap)
        end
    end

    ### save the leafs information to a CSV file
    partition = UEG.partition(maxOrder)
    leafStates = FeynmanDiagram.leafstates(leaf_maps, labelProd)
    len = length(leafStates)
    for (ikey, key) in enumerate(partition)
        key_str = join(string.(key))
        df = DataFrame([leafStates[idx][ikey] for idx in 1:len], :auto)
        CSV.write(root_dir * "leafinfo_GV$key_str.csv", df)
    end

    ### save the loop basis to a CSV file for the given maximum order
    df = DataFrame(labelProd.labels[end], :auto)
    CSV.write(root_dir * "loopBasis_GVmaxOrder$maxOrder.csv", df)
end

@inline function phase(varT, extT, l, β)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return exp(1im * π * (2l + 1) / β * (tout - tin))
end

#include("sigma_generic.jl")
include("sigmaKW.jl")
include("sigmaCuba.jl")
include("sigmaVegas.jl")
include("sigmaGV.jl")

function MC(para; kgrid=[para.kF,], ngrid=[-1, 0, 1], neval=1e6, reweight_goal=nothing,
    spinPolarPara::Float64=0.0, # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
    filename::Union{String,Nothing}=nothing, partition=UEG.partition(para.order), diagtype=:Parquet,
    isLayered2D=false # whether to use the screened Coulomb interaction in 2D or not 
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
        if spinPolarPara == 0.0
            diagram = Sigma.diagramGV(para, partition)
            sigma, result = Sigma.GV_unpolarized(para, diagram;
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
        diagram = Sigma.diagram(para, partition)
        sigma, result = Sigma.KW(para, diagram;
            neighbor=neighbor, reweight_goal=reweight_goal,
            kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread)
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

include("source_codeGV/func_sigmaGV_o1.jl")
include("source_codeGV/func_sigmaGV_o2.jl")
include("source_codeGV/func_sigmaGV_o3.jl")
include("source_codeGV/func_sigmaGV_o4.jl")
include("source_codeGV/func_sigmaGV_o5.jl")
const evalfunc_map = Dict(
    (1, 0, 0) => eval_graph100!,
    (1, 0, 1) => eval_graph101!,
    (1, 0, 2) => eval_graph102!,
    (1, 0, 3) => eval_graph103!,
    (1, 0, 4) => eval_graph104!,
    (1, 1, 0) => eval_graph110!,
    (1, 1, 1) => eval_graph111!,
    (1, 1, 2) => eval_graph112!,
    (1, 1, 3) => eval_graph113!,
    (1, 2, 0) => eval_graph120!,
    (1, 2, 1) => eval_graph121!,
    (1, 2, 2) => eval_graph122!,
    (1, 3, 0) => eval_graph130!,
    (1, 3, 1) => eval_graph131!,
    (1, 4, 0) => eval_graph140!,
    (2, 0, 0) => eval_graph200!,
    (2, 0, 1) => eval_graph201!,
    (2, 0, 2) => eval_graph202!,
    (2, 0, 3) => eval_graph203!,
    (2, 1, 0) => eval_graph210!,
    (2, 1, 1) => eval_graph211!,
    (2, 1, 2) => eval_graph212!,
    (2, 2, 0) => eval_graph220!,
    (2, 2, 1) => eval_graph221!,
    (2, 3, 0) => eval_graph230!,
    (3, 0, 0) => eval_graph300!,
    (3, 0, 1) => eval_graph301!,
    (3, 0, 2) => eval_graph302!,
    (3, 1, 0) => eval_graph310!,
    (3, 1, 1) => eval_graph311!,
    (3, 2, 0) => eval_graph320!,
    (4, 0, 0) => eval_graph400!,
    (4, 0, 1) => eval_graph401!,
    (4, 1, 0) => eval_graph410!,
    (5, 0, 0) => eval_graph500!,
    )

function set_evalfunc_map(maxorder::Int=5; root_dir=joinpath(@__DIR__, "source_codeGV/"))
    if maxorder == 6
        include(root_dir * "func_sigmaGV_o6.jl")
        evalfunc_map[(1,0,5)] = eval_graph105!
        evalfunc_map[(1,1,4)] = eval_graph114!
        evalfunc_map[(1,2,3)] = eval_graph123!
        evalfunc_map[(1,3,2)] = eval_graph132!
        evalfunc_map[(1,4,1)] = eval_graph141!
        evalfunc_map[(1,5,0)] = eval_graph150!
        evalfunc_map[(2,0,4)] = eval_graph204!
        evalfunc_map[(2,1,3)] = eval_graph213!
        evalfunc_map[(2,2,2)] = eval_graph222!
        evalfunc_map[(2,3,1)] = eval_graph231!
        evalfunc_map[(2,4,0)] = eval_graph240!
        evalfunc_map[(3,0,3)] = eval_graph303!
        evalfunc_map[(3,1,2)] = eval_graph312!
        evalfunc_map[(3,2,1)] = eval_graph321!
        evalfunc_map[(3,3,0)] = eval_graph330!
        evalfunc_map[(4,0,2)] = eval_graph402!
        evalfunc_map[(4,1,1)] = eval_graph411!
        evalfunc_map[(4,2,0)] = eval_graph420!
        evalfunc_map[(5,0,1)] = eval_graph501!
        evalfunc_map[(5,1,0)] = eval_graph510!
        evalfunc_map[(6,0,0)] = eval_graph600!
    end
end

end
