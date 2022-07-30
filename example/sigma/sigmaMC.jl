using ElectronLiquid
using JLD2

rs = [5.0,]
mass2 = [0.01,]
Fs = [0.0,]
beta = [25.0,]
order = [2,]
neval = 1e6

# mission = :Z
# mission = :SigmaK
mission = ARGS[1]
println("mission: ", mission)
# exit(0)

for _rs in rs
    for _mass2 in mass2
        for _F in Fs
            for _beta in beta
                for _order in order
                    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2)
                    kF = para.kF

                    if mission == "Z"
                        ######### calcualte Z factor ######################
                        kgrid = [kF,]
                        ngrid = [0, 1]
                    elseif mission == "K"
                        ######### calculate K dependence #####################
                        Nk, korder = 4, 4
                        minK = 0.2kF
                        kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 3kF], [kF,], Nk, minK, korder).grid
                        ngrid = [0,]
                    else
                        error("unknown mission")
                    end

                    sigma = ElectronLiquid.Sigma.sigmaKW(para, kgrid=kgrid, ngrid=ngrid, neval=neval)

                    jldopen("data_$(mission).jld2", "a+") do f
                        key = "$(UEG.short(para))"
                        if haskey(f, key)
                            @warn("replacing existing data for $key")
                            delete!(f, key)
                        end
                        f[key] = (para, ngrid, kgrid, sigma)
                    end
                end
            end
        end
    end
end

