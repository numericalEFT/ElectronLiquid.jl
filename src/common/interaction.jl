import ..ElectronGas: Interaction as Inter
import ..ElectronGas: Polarization
using ..ElectronGas: Parameter
using ..Lehmann

export Coulombinstant, KOinstant, KOstatic, interactionDynamic, interactionStatic, counterR

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

"""
    function KO_W(q, n, p::ParaMC)

KO interaction in momentum q and the Matsubara frequency index n
Assume
```math
r_q = v_q + f
```
then the KO interaction is
```math
Rq = r_q / (1 - r_q Π0) - f
```
"""
function KO_W(q, n::Integer, para::ParaMC; Pi=para.spin * Polarization.Polarization0_ZeroTemp(q, n, para.basic) * para.massratio)
    if abs(q) < 1e-6
        q = 1e-6
    end
    invKOinstant = 1.0 / KOinstant(q, para)
    Rs = 1.0 ./ (invKOinstant - Pi) - para.fs
    return Rs
end

function KOdynamic_T(para::ParaMC)
    return KOdynamic_T(para.basic, para.qgrid, para.τgrid, para.mass2, para.massratio, para.fs, para.fa)
end

function KOdynamic_T(basic::Parameter.Para, qgrid, τgrid, mass2=1.0e-6, massratio=1.0, fp=0.0, fm=0.0)
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

"""
    function counterKO_W(para::ParaMC; qgrid=para.qgrid, ngrid=[0,], order=para.order, proper=false)
    
    calculate counter-terms of the KO interaction

Assume
```math
r_q = v_q + f
```
and
```math
Rq = r_q / (1 - r_q Π0) - f
```
Then, the counter-term is given by a power expansion of the form
```math
(Rq ξ + f(ξ))/(1+(Rq ξ + f(ξ))Π0) - f(ξ)
```
where f(ξ) = f1 ξ + f2 ξ^2 + ...

Therefore, the counter-term is given by
Order 1 (ξ^2)
```math
(Rq+f1)^2 Π0
```
Order 2 (ξ^3)
```math
(Rq+f1)^3 Π0^2 + (R_q+f1)f2 Π_0
```

"""
function counterKO_W(para::ParaMC; qgrid=para.qgrid, ngrid=[0,], order=para.order, proper=false, bubble=true)
    Nq, Nw = length(qgrid), length(ngrid)
    cRs1 = zeros(Float64, (Nq, Nw))
    for (ni, n) in enumerate(ngrid)
        for (qi, q) in enumerate(qgrid)
            Pi = para.spin * Polarization.Polarization0_ZeroTemp(q, n, para.basic) * para.massratio
            if proper == false
                Rs = KO_W(q, n, para; Pi=Pi)
            else
                Rs = 0.0
            end
            # Rs will be zero for the proper counter-term
            if order == 1
                cRs1[qi, ni] = -(Rs + para.fs)^2 * Pi
            elseif order == 2
                cRs1[qi, ni] = (Rs + para.fs)^3 * Pi^2
            else
                error("not implemented!")
            end

            if bubble == false
                cRs1[qi, ni] += (-Rs)^(order + 1) * Pi^order
            end
        end
    end
    return cRs1
end

function counterKO_T(para::ParaMC; qgrid=para.qgrid, τgrid=para.τgrid, order=para.order, proper=false, bubble=true)
    dlr = DLRGrid(Euv=10 * para.EF, β=para.β, rtol=1e-10, isFermi=false, symmetry=:ph) # effective interaction is a correlation function of the form <O(τ)O(0)>
    cRs1 = counterKO_W(para; qgrid=qgrid, ngrid=dlr.n, order=order, proper=proper, bubble=bubble)
    cRs1 = matfreq2tau(dlr, cRs1, τgrid.grid, axis=2)
    return real.(cRs1)
end

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
    qgrid, τgrid = p.qgrid, p.τgrid
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

"""
    function interactionDynamic(p::ParaMC, qd, τIn, τOut)

Dynamic part of the interaction.

Assume 
```math
r_q = v_q + f
```
Then, the dynamic interaction is given by

```math
δR_q = (r_q)²Π₀/(1-r_q Π₀)
```

where Π₀ is the polarization of free electrons.
"""
function interactionDynamic(p::ParaMC, qd, τIn, τOut)
    # @unpack qgrid, τgrid = p.qgrid, p.τgrid
    @unpack qgrid, τgrid = p
    kF, maxK, dW0 = p.kF, p.maxK, p.dW0

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
    # println(qgrid)
    # println(τgrid)
    # exit(0)
    return vd * linear2D(dW0, qgrid, τgrid, qd, dτ) # dynamic interaction, don't forget the singular factor vq
end

"""
    function interactionStatic(p::ParaMC, qd, τIn, τOut)
    
instant part of the renormalized interaction


Assume 
```math
r_q = v_q + f
```
Then, the instant interaction is given by 
```math
v_q = r_q - f
```

The current implementation involves one auxiliary time variable τOut for better sign cancellation.

To show the net result is v_q, one may perform a τOut integration explicitly, then
```math
kostatic = r_q / (1-r_q Π₀) - f
```
where Π₀ is the polarization of free electrons.
and,
```math
dynamic =  (r_q)²Π₀/(1-r_q Π₀)
```
so that,
```math
kostatic - dynamic = v_q
```
"""
function interactionStatic(p::ParaMC, qd, τIn, τOut)
    kF, maxK, β = p.kF, p.maxK, p.β

    if qd > maxK
        return 0.0
    end
    if qd <= 1e-6 * kF
        qd = 1e-6 * kF
    end
    # if there is no dynamic interactoin
    # return KOinstant(qd)

    # one must divide by beta because there is an auxiliary time variable for each interaction
    # return KOinstant(qd, p) / β

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
