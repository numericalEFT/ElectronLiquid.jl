@testset "ParaMC" begin
    # Test that ParaMC can be converted to a string and back
    para1 = UEG.ParaMC(rs=2.0, beta=25.0, order=4, mass2=0.5, isDynamic=false, additional=[1, 2, 3])
    s = UEG.short(para1)
    @test s == "Fa_-0.0_Fs_-0.0_beta_25.0_dim_3_isDynamic_false_isFock_false_mass2_0.5_massratio_1.0_order_4_rs_2.0_spin_2"
    para1p = UEG.ParaMC(s)
    @test isequal(para1p, para1)
    @test para1p == para1
    # Additional data is lost after conversion and ignored for purposes of equality
    @test para1p.additional == []
    @test para1p.additional != para1.additional

    # Test equality before/after initialization
    para2 = UEG.ParaMC(rs=5.0, beta=25.0, Fs=-1.0, order=1, mass2=1e-5, isDynamic=true)
    para3 = UEG.ParaMC(rs=5.0, beta=25.0, Fs=-1.0, order=1, mass2=1e-5, isDynamic=true)
    para4 = UEG.ParaMC(rs=5.0, beta=25.0, Fs=-1.5, order=1, mass2=1e-5, isDynamic=true)
    # NOTE: we may have para2.dW0 != para3.dW0 before initialization
    @test para2 != para1
    @test para2 == para3
    @test para2 != para4
    # Initialize dynamic ParaMC objects
    UEG.MCinitialize!(para2)
    UEG.MCinitialize!(para3)
    UEG.MCinitialize!(para4)
    # Test equality after initialization (include initialized fields)
    @test para2 != para1
    @test para2 == para3
    @test para2 != para4
end
