@testset "RG_treelevel_interaction" begin
    include("../example/RG/gamma4_treelevel.jl")
    para = UEG.ParaMC(rs=4.0, beta=100.0, order=1, mass2=1e-5, isDynamic=true, isFock=false, Fs=-0.0, Fa=-0.0)

    Λgrid = CompositeGrid.LogDensedGrid(:gauss, [1.0 * para.kF, 100 * para.kF], [para.kF,], 8, 0.01 * para.kF, 8)

    fs, us, dfs, dus = gamma4_treelevel_KO(para, Λgrid)

    _fs = -Interp.integrate1D(dfs, Λgrid, [Λgrid[1], Λgrid[end]])

    @test abs(_fs - fs[1]) < 1e-3

    paras = [UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs[li], Fa=0.0, order=para.order,
        mass2=para.mass2, isDynamic=true, isFock=false) for li in eachindex(Λgrid)]

    # d R_ex/df on FS
    dR = ∂Rs_∂fs_exchange(paras, [para.kF for l in Λgrid]; ct=false) / para.NF

    # R(f=0) - R(f=f_0)
    diff = Interp.integrate1D(dR .* dfs, Λgrid, [Λgrid[1], Λgrid[end]])

    wp, wm, angle = Ver4.exchange_interaction(para, para.kF, para.kF; ct=false)
    W0 = Ver4.Legrendre(0, wp, angle)

    wp, wm, angle = Ver4.exchange_interaction(paras[1], para.kF, para.kF; ct=false)
    R_f = Ver4.Legrendre(0, wp, angle)

    println(" W0: ", W0, " R_f: ", R_f, " diff: ", diff)
    @test abs(diff - (W0 - R_f)) < 1e-3
end