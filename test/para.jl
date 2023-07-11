@testset "ParaMC" begin
    # Test that ParaMC can be converted to a string and back
    para = UEG.ParaMC(rs=2.0, beta=25.0, order=4, mass2=0.5, isDynamic=false)
    s = UEG.short(para)
    @test s == "Fa_-0.0_Fs_-0.0_beta_25.0_dim_3_isDynamic_false_isFock_false_mass2_0.5_massratio_1.0_order_4_rs_2.0_spin_2"
    para_from_s = UEG.ParaMC(s)
    @test isequal(para_from_s, para)
    @test para_from_s == para

    # Test equality before/after initialization for dynamic case
    mass2 = 1e-5
    Fs = -1.0
    dF = 0.0001
    para1 = UEG.ParaMC(rs=5.0, beta=25.0, Fs=Fs, order=1, mass2=mass2, isDynamic=true)
    para2 = UEG.ParaMC(rs=5.0, beta=25.0, Fs=Fs + dF, order=1, mass2=mass2, isDynamic=true)
    # Parameters are equal before initialization (ignore uninitialized fields)
    @test para1 == para2
    UEG.MCinitialize!(para1)
    UEG.MCinitialize!(para2)
    # Parameters are equal after initialization (include uninitialized fields)
    @test para1 == para2
end
