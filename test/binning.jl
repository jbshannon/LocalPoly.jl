@testset "Binning" begin
    # 1-dimensional
    y = 2
    grid = (g = 0.0:1.0:1.0, c = zeros(2), d = zeros(2))
    for x in 0.2:0.1:0.8
        LocalPoly.linear_binning!(grid, [x], [y])
        @test grid.c[1] ≈ 1-x
        @test grid.c[2] ≈ x
        @test grid.d ≈ grid.c*y
        fill!(grid.c, 0); fill!(grid.d, 0);
    end

    # 2-dimensional
    y = 2
    grid = (g = (0.0:1.0:1.0, 0.0:1.0:1.0), c = zeros(2, 2), d = zeros(2, 2))
    for x in Iterators.product(0.2:0.2:0.8, 0.2:0.2:0.8)
        LocalPoly.linear_binning!(grid, hcat(x...), [y])
        @test grid.c[1, 1] ≈ (1 - x[1]) * (1 - x[2])
        @test grid.c[1, 2] ≈ (1 - x[1]) * x[2]
        @test grid.c[2, 1] ≈ x[1] * (1 - x[2])
        @test grid.c[2, 2] ≈ x[1] * x[2]
        @test grid.d ≈ grid.c*y
        fill!(grid.c, 0); fill!(grid.d, 0);
    end

    # 3-dimensional
    y = 2
    g = ntuple(i -> 0.0:1.0, 3)
    M = length.(g)
    grid = (g = g, c = zeros(M), d = zeros(M))
    x = [0.3 0.4 0.5]
    LocalPoly.linear_binning!(grid, x, [y])
    @test grid.c[1, 1, 1] ≈ (1 - x[1]) * (1 - x[2]) * (1 - x[3])
    @test grid.c[2, 1, 1] ≈       x[1] * (1 - x[2]) * (1 - x[3])
    @test grid.c[1, 2, 1] ≈ (1 - x[1]) *       x[2] * (1 - x[3])
    @test grid.c[1, 1, 2] ≈ (1 - x[1]) * (1 - x[2]) * x[3]
    @test grid.c[2, 2, 1] ≈       x[1] *       x[2] * (1 - x[3])
    @test grid.c[2, 1, 2] ≈       x[1] * (1 - x[2]) * x[3]
    @test grid.c[1, 2, 2] ≈ (1 - x[1]) *       x[2] * x[3]
    @test grid.c[2, 2, 2] ≈       x[1] *       x[2] * x[3]
end
