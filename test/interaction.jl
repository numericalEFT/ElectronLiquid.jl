@testset "Interaction" begin
    mass2 = 1e-5
    Fs = -1.0
    dF = 0.0001
    ####################### STATIC #########################
    para1 = UEG.ParaMC(rs=5.0, beta=25.0, Fs=Fs, order=1, mass2=mass2, isDynamic=true)
    UEG.MCinitialize!(para1)

    para2 = UEG.ParaMC(rs=5.0, beta=25.0, Fs=Fs + dF, order=1, mass2=mass2, isDynamic=true)
    UEG.MCinitialize!(para2)

    # dW0_df = (para2.dW0 - para1.dW0) / (dF / para1.NF)
    # dW0 = (para1.dW0 + para2.dW0) / 2
    # avg = (para1.dW0_f + para2.dW0_f) / 2
    # for (qi, q) in enumerate(para1.qgrid.grid)
    #     for (τi, τ) in enumerate(para1.τgrid.grid)
    #         dW0_df[qi, τi] = dW0_df[qi, τi] * UEG.KOinstant(q, para1) + dW0[qi, τi]
    #     end
    # end
    # diff = dW0_df - avg
    # @test maximum(abs.(diff)) < 1e-5
    # println(dW0_df[1:10])
    # println(avg[1:10])
    # println(diff[1:10])
end