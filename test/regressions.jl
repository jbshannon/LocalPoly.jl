@testset verbose=true "Regressions" begin
    x = 2π * rand(100_000)
    y = sin.(x) + randn(size(x))/10
    v = range(minimum(x), maximum(x), length=1000)

    # Test that function value is well-approximated
    β̂ = lpreg(y, x, v; nbins=1000)
    I = findall(β -> !isnan(β[1]), β̂)
    sqerr = map(i -> abs2(β̂[i][1] - sin(v[i])), I)
    @test mean(sqerr) < 1.0 # what would a good threshold be?

    # Test that first derivative is well-approximated
    β̂ = lpreg(y, x, v; degree=2, nbins=1000)
    I = findall(β -> !isnan(β[2]), β̂)
    sqerr = map(i -> abs2(β̂[i][2] - cos(v[i])), I)
    @test mean(sqerr) < 1.0 # what would a good threshold be?
end
