@testset "Regressions" begin
    x = 2π * rand(100_000)
    y = sin.(x) + randn(size(x))/10
    v = range(minimum(x), maximum(x), length=1000)

    # Test that function value is well-approximated
    β̂ = lpreg(x, y, v; nbins=1000)
    MSE = mean(map((b, x) -> abs2(b[1] - sin(x)), β̂, v))
    @test MSE < 1.0 # what would a good threshold be?

    # Test that first derivative is well-approximated
    β̂ = lpreg(x, y, v; ν=1, degree=2, nbins=1000)
    MSE = mean(map((b, x) -> abs2(b[2] - cos(x)), β̂, v))
    @test MSE < 1.0 # what would a good threshold be?
end
