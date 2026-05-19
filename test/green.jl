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
