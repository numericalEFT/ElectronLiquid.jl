function integrandKW(idx, vars, config)
    varK, varT, varN, ExtKidx = vars
    para, kgrid, ngrid, MaxLoopNum, extT_labels = config.userdata[1:5]
    leafval, leafType, leafOrders, leafτ_i, leafτ_o, leafMomIdx = config.userdata[6]
    momLoopPool, root = config.userdata[7:8]
    graphfuncs! = config.userdata[9][idx]
    isLayered2D = config.userdata[end]
    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0

    extidx = ExtKidx[1]
    varK.data[1, 1] = kgrid[extidx]
    FrontEnds.update(momLoopPool, varK.data[:, 1:MaxLoopNum])
    for (i, lftype) in enumerate(leafType[idx])
        if lftype == 0
            continue
            # elseif isodd(lftype) #fermionic 
        elseif lftype == 1 #fermionic 
            τ = varT[leafτ_o[idx][i]] - varT[leafτ_i[idx][i]]
            kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
            ϵ = dot(kq, kq) / (2me) - μ
            order = leafOrders[idx][i][1]
            leafval[idx][i] = Propagator.green_derive(τ, ϵ, β, order)
        elseif lftype == 2 #bosonic 
            kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
            # order = lftype / 2 - 1
            order = leafOrders[idx][i][2]
            if dim == 3
                invK = 1.0 / (dot(kq, kq) + λ)
                leafval[idx][i] = e0^2 / ϵ0 * invK * (λ * invK)^order
            elseif dim == 2
                if isLayered2D == false
                    invK = 1.0 / (sqrt(dot(kq, kq)) + λ)
                    leafval[idx][i] = e0^2 / 2ϵ0 * invK * (λ * invK)^order
                else
                    if order == 0
                        q = sqrt(dot(kq, kq) + 1e-16)
                        invK = 1.0 / q
                        leafval[idx][i] = e0^2 / 2ϵ0 * invK * tanh(λ * q)
                    else
                        leafval[idx][i] = 0.0 # no high-order counterterms
                    end
                end
            else
                error("not implemented!")
            end
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    graphfuncs!(root, leafval[idx])  # allocations due to run-time variable `idx`

    n = ngrid[varN[1]]
    weight = sum(root[i] * phase(varT, extT, n, β) for (i, extT) in enumerate(extT_labels[idx]))

    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(dim * loopNum)
    return weight * factor
end

function integrandKW_Clib(idx, vars, config)
    varK, varT, varN, ExtKidx = vars
    para, kgrid, ngrid, MaxLoopNum, extT_labels = config.userdata[1:5]
    leafstates, leafval = config.userdata[6][idx], config.userdata[7][idx]
    momLoopPool, root = config.userdata[8:9]
    isLayered2D = config.userdata[10]
    partition = config.userdata[11]

    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0
    extidx = ExtKidx[1]
    varK.data[1, 1] = kgrid[extidx]

    FrontEnds.update(momLoopPool, varK.data[:, 1:MaxLoopNum])
    for (i, lfstat) in enumerate(leafstates)
        lftype, lforders, leafτ_i, leafτ_o, leafMomIdx = lfstat.type, lfstat.orders, lfstat.inTau_idx, lfstat.outTau_idx, lfstat.loop_idx
        if lftype == 0
            continue
            # elseif isodd(lftype) #fermionic 
        elseif lftype == 1 #fermionic 
            τ = varT[leafτ_o] - varT[leafτ_i]
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            ϵ = dot(kq, kq) / (2me) - μ
            # order = (lftype - 1) / 2
            order = lforders[1]
            leafval[i] = Propagator.green_derive(τ, ϵ, β, order)
        elseif lftype == 2 #bosonic
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            # order = lftype / 2 - 1
            order = lforders[2]
            if dim == 3
                invK = 1.0 / (dot(kq, kq) + λ)
                leafval[i] = e0^2 / ϵ0 * invK * (λ * invK)^order
            elseif dim == 2
                if isLayered2D == false
                    invK = 1.0 / (sqrt(dot(kq, kq)) + λ)
                    leafval[i] = e0^2 / 2ϵ0 * invK * (λ * invK)^order
                else
                    if order == 0
                        q = sqrt(dot(kq, kq) + 1e-16)
                        invK = 1.0 / q
                        leafval[i] = e0^2 / 2ϵ0 * invK * tanh(λ * q)
                    else
                        leafval[i] = 0.0 # no high-order counterterms
                    end
                end
            else
                error("not implemented!")
            end
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    group = partition[idx]
    evalfuncParquetAD_sigma_map[group](root, leafval)

    n = ngrid[varN[1]]
    weight = sum(root[i] * phase(varT, extT, n, β) for (i, extT) in enumerate(extT_labels[idx]))
    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(dim * loopNum)
    return weight * factor
end