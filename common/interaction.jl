import ElectronGas: Interaction as Inter
import ElectronGas: Polarization
using ElectronGas: Parameter
using Lehmann

"""
    function lindhard(x)

    Dimensionless Linhard function so that lindhard(0) = 1 for both 2D and 3D

# Arguments

- x : 2k/k_F
"""
@inline function lindhard(x, dim)
    if dim == 3
        if (abs(x) < 1.0e-4)
            return 1.0
        elseif (abs(x - 1.0) < 1.0e-4)
            return 0.5
        else
            return 0.5 - (x^2 - 1) / 4.0 / x * log(abs((1 + x) / (1 - x)))
        end
    elseif dim == 2
        if (abs(x) <= 1.0)
            return 1.0
        else
            return 1 - sqrt(1 - 1 / x^2)
        end
    else
        error("not implemented!")
    end
end

@inline function Coulombinstant(q, p::ParaMC)
    return Coulombinstant(q, p.basic, p.mass2)
end

@inline function Coulombinstant(q, basic, mass2=1e-6)
    e0, dim = basic.e0, basic.dim
    if abs(q) < 1e-6
        q = 1e-6
    end
    if dim == 3
        return 4π * e0^2 / (q^2 + mass2)
    elseif dim == 2
        return 2π * e0^2 / sqrt(q^2 + mass2)
    else
        error("not implemented!")
    end
end

@inline function KOinstant(q, p::ParaMC)
    return KOinstant(q, p.basic, p.mass2, p.massratio, p.fs, p.fa)
end

@inline function KOinstant(q, basic, mass2=1e-6, massratio=1.0, fp=0.0, fm=0.0)
    e0, dim = basic.e0, basic.dim
    if abs(q) < 1e-6
        q = 1e-6
    end
    if dim == 3
        return 4π * e0^2 / (q^2 + mass2) + fp
    elseif dim == 2
        return 2π * e0^2 / sqrt(q^2 + mass2) + fp
    else
        error("not implemented!")
    end
end

@inline function KOstatic(q, p::ParaMC)
    return KOstatic(q, p.basic, p.mass2, p.massratio, p.fs, p.fa)
end

@inline function KOstatic(q, basic, mass2=1e-6, massratio=1.0, fp=0.0, fm=0.0)
    e0, dim = basic.e0, basic.dim
    kF, NF = basic.kF, basic.NF

    Pi = -lindhard(q / 2.0 / kF, dim) * NF * massratio
    if dim == 3
        vd = (4π * e0^2 + fp * (q^2 + mass2)) / ((1 - fp * Pi) * (q^2 + mass2) - 4π * e0^2 * Pi) - fp
        return vd
    elseif dim == 2
        qm = sqrt(q^2 + mass2)
        vd = (2π * e0^2 + fp * qm) / ((1 - fp * Pi) * qm - 2π * e0^2 * Pi) - fp
        return vd
    else
        error("not implemented!")
    end
end

function KO(basic::Parameter.Para, qgrid, τgrid, mass2=1.0e-6, massratio=1.0, fp=0.0, fm=0.0)
    # para = Parameter.rydbergUnit(1.0 / beta, rs, dim, Λs=mass2)
    dim, e0 = basic.dim, basic.e0
    dlr = DLRGrid(Euv=10 * basic.EF, β=basic.β, rtol=1e-10, isFermi=false, symmetry=:ph) # effective interaction is a correlation function of the form <O(τ)O(0)>
    Nq, Nτ = length(qgrid), length(τgrid)
    Rs = zeros(Float64, (Nq, dlr.size)) # Matsubara grid is the optimized sparse DLR grid 
    Ra = zeros(Float64, (Nq, dlr.size)) # Matsubara grid is the optimized sparse DLR grid 
    Pi = zeros(Float64, (Nq, dlr.size)) # Matsubara grid is the optimized sparse DLR grid 
    for (ni, n) in enumerate(dlr.n)
        for (qi, q) in enumerate(qgrid)
            invKOinstant = 1.0 / KOinstant(q, basic, mass2, massratio, fp, fm)
            # Rs = (vq+f)Π0/(1-(vq+f)Π0)
            Pi[qi, ni] = basic.spin * Polarization.Polarization0_ZeroTemp(q, n, basic) * massratio
            Rs[qi, ni] = Pi[qi, ni] / (invKOinstant - Pi[qi, ni])
            # if dim == 3
            #     a = Inter.KO(q, n, para, landaufunc=Inter.landauParameterConst,
            #         Fs=-Fs, Fa=-Fa, massratio=massratio, regular=true)[1]
            #     @assert abs(Rs[qi, ni] - a) < 1e-10 "$(Rs[qi, ni]) vs $a with diff = $(Rs[qi, ni]-a)"
            # end
        end
    end
    for (qi, q) in enumerate(qgrid)
        # turn on this to check consistencey between two static KO interactions
        w = KOinstant(q, basic, mass2, massratio, fp, fm) * Rs[qi, 1] + Coulombinstant(q, basic, mass2)
        kostatic = KOstatic(q, basic, mass2, massratio, fp, fm)
        @assert abs(w - kostatic) < 1e-4 "$q  ==> $w != $(kostatic)"
        # println("$(q/kF)   $(w*NF)")
        staticPi = -basic.NF * massratio * lindhard(q / 2 / basic.kF, basic.dim)
        @assert abs(Pi[qi, 1] - staticPi) < 1e-8 "$(Pi[qi, 1]) vs $staticPi"
    end
    # exit(0)
    # println(Rs[:, 1])
    Rs = matfreq2tau(dlr, Rs, τgrid.grid, axis=2)
    # for (qi, q) in enumerate(qgrid)
    #     println("$(q/kF)   $(Rs[qi, 1])")
    # end
    return real.(Rs)
end

# const dW0 = KO(qgrid, τgrid)

function counterKO(basic::Parameter.Para, qgrid, τgrid, order, mass2=1.0e-6, massratio=1.0, fp=0.0, fm=0.0)
    # para = Parameter.rydbergUnit(1.0 / beta, rs, dim, Λs=mass2)
    dim, e0 = basic.dim, basic.e0
    dlr = DLRGrid(Euv=10 * basic.EF, β=basic.β, rtol=1e-10, isFermi=false, symmetry=:ph) # effective interaction is a correlation function of the form <O(τ)O(0)>
    Nq, Nτ = length(qgrid), length(τgrid)
    Rs = zeros(Float64, (Nq, dlr.size)) # Matsubara grid is the optimized sparse DLR grid 
    Ra = zeros(Float64, (Nq, dlr.size)) # Matsubara grid is the optimized sparse DLR grid 
    Pi = zeros(Float64, (Nq, dlr.size)) # Matsubara grid is the optimized sparse DLR grid 
    cRs1 = zeros(Float64, (Nq, dlr.size))
    for (ni, n) in enumerate(dlr.n)
        for (qi, q) in enumerate(qgrid)
            invKOinstant = 1.0 / KOinstant(q, basic, mass2, massratio, fp, fm)
            # Rs = (vq+f)Π0/(1-(vq+f)Π0)
            Pi[qi, ni] = basic.spin * Polarization.Polarization0_ZeroTemp(q, n, basic) * massratio
            Rs[qi, ni] = Pi[qi, ni] / (invKOinstant - Pi[qi, ni])
            cRs1[qi, ni] = (-Rs[qi, ni])^order / (invKOinstant - Pi[qi, ni])
        end
    end
    cRs1 = matfreq2tau(dlr, cRs1, τgrid.grid, axis=2)
    return real.(cRs1)
end

# const cRs1 = counterKO(qgrid, τgrid, 1)
# const cRs2 = counterKO(qgrid, τgrid, 2)

"""
   linear2D(data, xgrid, ygrid, x, y) 

linear interpolation of data(x, y)

#Arguments:
- xgrid: one-dimensional grid of x
- ygrid: one-dimensional grid of y
- data: two-dimensional array of data
- x: x
- y: y
"""
@inline function linear2D(data, xgrid, ygrid, x, y)

    xarray, yarray = xgrid.grid, ygrid.grid

    xi0, xi1, yi0, yi1 = 0, 0, 0, 0
    if (x <= xarray[firstindex(xgrid)])
        xi0 = 1
        xi1 = 2
    elseif (x >= xarray[lastindex(xgrid)])
        xi0 = lastindex(xgrid) - 1
        xi1 = xi0 + 1
    else
        xi0 = floor(xgrid, x)
        xi1 = xi0 + 1
    end

    if (y <= yarray[firstindex(ygrid)])
        yi0 = 1
        yi1 = 2
    elseif (y >= yarray[lastindex(ygrid)])
        yi0 = lastindex(ygrid) - 1
        yi1 = yi0 + 1
    else
        yi0 = floor(ygrid, y)
        yi1 = yi0 + 1
    end

    dx0, dx1 = x - xarray[xi0], xarray[xi1] - x
    dy0, dy1 = y - yarray[yi0], yarray[yi1] - y

    g0 = data[xi0, yi0] * dx1 + data[xi1, yi0] * dx0
    g1 = data[xi0, yi1] * dx1 + data[xi1, yi1] * dx0

    gx = (g0 * dy1 + g1 * dy0) / (dx0 + dx1) / (dy0 + dy1)
    return gx
end

function counterR(p::ParaMC, qd, τIn, τOut, order)
    kF, maxK = p.kF, p.maxK
    qgrid, τgrid = q.qgrid, q.τgrid
    cRs = p.cRs

    if qd > maxK
        return 0.0
    end

    dτ = abs(τOut - τIn)

    # if qd <= qgrid.grid[1]
    # the current interpolation vanishes at q=0, which needs to be corrected!
    if qd <= 1e-6 * kF
        # q = qgrid.grid[1] + 1.0e-6
        qd = 1e-6 * kF
    end

    if order <= p.order
        return linear2D(cRs[order], qgrid, τgrid, qd, dτ)
    else
        error("not implemented!")
    end
end

function interactionDynamic(p::ParaMC, qd, τIn, τOut)
    dim, e0, kF = p.dim, p.e0, p.kF
    qgrid, τgrid = q.qgrid, q.τgrid
    dW0 = p.dW0

    if qd > maxK
        return 0.0
    end

    dτ = abs(τOut - τIn)

    # if qd <= qgrid.grid[1]
    # the current interpolation vanishes at q=0, which needs to be corrected!
    if qd <= 1e-6 * kF
        # q = qgrid.grid[1] + 1.0e-6
        qd = 1e-6 * kF
    end

    vd = KOinstant(qd, p)
    return vd * linear2D(dW0, qgrid, τgrid, qd, dτ) # dynamic interaction, don't forget the singular factor vq
end

function interactionStatic(p::ParaMC, qd, τIn, τOut)
    kF, maxK, β = p.kF, p.maxK, p.β
    dim, e0, NF = p.dim, p.e0, p.NF

    if qd > maxK
        return 0.0
    end
    if qd <= 1e-6 * kF
        qd = 1e-6 * kF
    end
    # if there is no dynamic interactoin
    # return KOinstant(qd)

    # one must divide by beta because there is an auxiliary time variable for each interaction
    # return KOinstant(qd) / β

    # introduce a fake tau variable to alleviate sign cancellation between the static and the dynamic interactions
    # if qd > 50 * kF
    #     println("$τIn, $τOut")
    #     println("$(KOstatic(qd) / β), $(interactionDynamic(qd, τIn, τOut)), $(fp / β)")
    #     exit(0)
    # end
    kostatic = KOstatic(qd, p)
    return kostatic / β - interactionDynamic(p, qd, τIn, τOut)
end

# const qgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 6 * kF], [0.0, 2kF], 16, 0.01 * kF, 8)
# const τgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 16, β * 1e-4, 8)
# vqinv = [(q^2 + mass2) / (4π * e0^2) for q in qgrid.grid]
# const dW0 = TwoPoint.dWRPA(vqinv, qgrid.grid, τgrid.grid, dim, EF, kF, β, spin, me) # dynamic part of the effective interaction
