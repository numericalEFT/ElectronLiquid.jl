@testset "Green Parquet diagrams" begin
    para = ParaMC(; order=1, rs=1.0, beta=10.0, isDynamic=false)
    partition = [(0, 0, 0), (1, 0, 0)]

    diagram = Diagram.diagram_parquet_noresponse(
        :green,
        para,
        partition;
        filter=[FeynmanDiagram.Parquet.NoHartree],
        optimize_level=1,
    )

    @test diagram[1] == partition
    @test length(diagram) == 4
    @test all(length(extT) == length(diagram[3][p]) for (p, extT) in zip(diagram[1], diagram[4]))
    @test all(all(length(label) == 2 for label in extT) for extT in diagram[4])

    green_diagram = Green.diagram(
        para,
        partition;
        filter=[FeynmanDiagram.Parquet.NoHartree],
        optimize_level=1,
    )

    @test green_diagram[1] == partition
end

@testset "Green result output" begin
    para = ParaMC(; order=1, rs=1.0, beta=10.0, isDynamic=false)
    partition = [(1, 0, 0)]
    kgrid = [para.kF]
    tgrid = [para.beta / 2]
    ngrid = [0]

    kt_data = Dict(partition[1] => fill(measurement(1.25, 0.125), length(tgrid), length(kgrid)))
    kt_output = sprint() do io
        Green._print_data(io, para, partition, tgrid, kgrid, kt_data, :KT)
    end
    @test occursin("Group (1, 0, 0)", kt_output)
    @test occursin("t = $(tgrid[1])", kt_output)
    @test occursin("q/kF", kt_output)
    @test occursin("avg", kt_output)
    @test occursin("err", kt_output)
    @test occursin("  1.000000    1.250000 ±   0.125000", kt_output)

    kw_data = Dict(partition[1] => fill(
        complex(measurement(1.25, 0.125), measurement(-0.5, 0.05)),
        length(ngrid),
        length(kgrid),
    ))
    kw_output = sprint() do io
        Green._print_data(io, para, partition, ngrid, kgrid, kw_data, :KW)
    end
    @test occursin("n = 0", kw_output)
    @test occursin("real(avg)", kw_output)
    @test occursin("imag(avg)", kw_output)
    @test occursin("  1.000000    1.250000 ±   0.125000    -0.500000 ±   0.050000", kw_output)

    mktempdir() do dir
        filename = joinpath(dir, "green.jld2")
        Green._save_data(filename, para, tgrid, kgrid, kt_data)
        ElectronLiquid.Green.JLD2.jldopen(filename, "r") do f
            stored_grid, stored_kgrid, stored_data = f["$(UEG.short(para))"]
            @test stored_grid == tgrid
            @test stored_kgrid == kgrid
            @test stored_data == kt_data
        end
    end
end
