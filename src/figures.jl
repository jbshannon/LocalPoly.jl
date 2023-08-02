using CairoMakie

github_light = Theme(
    backgroundcolor = :transparent,
    textcolor = :gray50,
    Axis = (
        backgroundcolor = :transparent,
        xgridcolor = (:black, 0.07),
        ygridcolor = (:black, 0.07),
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
        xminorticksvisible = false,
        yminorticksvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        xlabelpadding = 3,
        ylabelpadding = 3
    ),
    Legend = (
        framevisible = false,
        padding = (0, 0, 0, 0),
    ),
    Axis3 = (
        xgridcolor = (:black, 0.07),
        ygridcolor = (:black, 0.07),
        zgridcolor = (:black, 0.07),
        xspinesvisible = false,
        yspinesvisible = false,
        zspinesvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        zticksvisible = false,
    ),
    Colorbar = (
        ticksvisible = false,
        spinewidth = 0,
        ticklabelpad = 5,
    )
)

github_dark = Theme(
    backgroundcolor = :transparent,
    textcolor = :gray45,
    linecolor = :gray60,
    Axis = (
        backgroundcolor = :transparent,
        xgridcolor = (:white, 0.09),
        ygridcolor = (:white, 0.09),
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
        xminorticksvisible = false,
        yminorticksvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        xlabelpadding = 3,
        ylabelpadding = 3
    ),
    Legend = (
        framevisible = false,
        padding = (0, 0, 0, 0),
    ),
    Axis3 = (
        xgridcolor = (:white, 0.09),
        ygridcolor = (:white, 0.09),
        zgridcolor = (:white, 0.09),
        xspinesvisible = false,
        yspinesvisible = false,
        zspinesvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        zticksvisible = false,
    ),
    Colorbar = (
        ticksvisible = false,
        spinewidth = 0,
        ticklabelpad = 5,
    )
)

##


##
figdir = joinpath(
    dirname(@__DIR__),
    "docs",
    "src",
    "images",
)

##
using LocalPoly, Random
Random.seed!(42)
x = 2π * rand(1000)
y = sin.(x) + randn(size(x))/4
grid = linear_binning(x, y; nbins=100)
v = gridnodes(grid)
β̂ = lpreg(grid)
ŷ = first.(β̂)
ci = confint(grid)

ν = 1
degree = ν+1
h = plugin_bandwidth(grid; ν)
β̂1 = lpreg(grid; degree, h)
ŷ′ = [b[ν+1] for b in β̂1]
ci1 = confint(grid; ν, p=degree)
##
using CairoMakie
fig, ax, sc = scatter!(x, y; markersize=5, label="Data")
lines!(ax, v, sin.(v); color=:darkgreen, label="True values")
lines!(ax, v, ŷ; color=:tomato, linewidth=3, label="Fitted values")
Legend(fig[2, 1], ax; orientation=:horizontal, framevisible=false)
current_figure()

##
function _fit()
    fig, ax, sc = scatter(x, y; color=Cycled(1), markersize=3, label="Data")
    lines!(ax, v, sin.(v); label="True values")
    lines!(ax, v, ŷ; linewidth=3, label="Fitted values")
    Legend(fig[2, 1], ax; orientation=:horizontal)
    return fig
end

lightfig = with_theme(_fit, github_light)
save(joinpath(figdir, "readme", "light", "fit.svg"), lightfig)
save(joinpath(figdir, "example", "fit.svg"), lightfig)

darkfig = with_theme(_fit, github_dark)
save(joinpath(figdir, "readme", "dark", "fit.svg"), darkfig)
@info "Saved fit"
##

function _ci()
    fig, ax, sc = scatter(x, y; color=Cycled(1), markersize=3, label="Data")
    lines!(ax, v, sin.(v); label="True values")
    lines!(ax, v, ŷ; linewidth=3, label="Fitted values")
    band!(ax, v, first.(ci), last.(ci); color=(:tomato, 0.3))
    Legend(fig[2, 1], ax; orientation=:horizontal)
    return fig
end

lightfig = with_theme(_ci, github_light)
save(joinpath(figdir, "readme", "light", "fit_ci.svg"), lightfig)
save(joinpath(figdir, "example", "fit_ci.svg"), lightfig)

darkfig = with_theme(_ci, github_dark)
save(joinpath(figdir, "readme", "dark", "fit_ci.svg"), darkfig)
@info "Saved fit with confidence interval"

##

function _derivative()
    fig, ax, sc = scatter(x, y; color=Cycled(1), markersize=3, label="Data")
    lines!(ax, v, cos.(v); label="True values")
    lines!(ax, v, ŷ′; linewidth=3, label="Fitted values")
    band!(ax, v, first.(ci1), last.(ci1); color=(:tomato, 0.3))
    Legend(fig[2, 1], ax; orientation=:horizontal)
    return fig
end

lightfig = with_theme(_ci, github_light)
save(joinpath(figdir, "readme", "light", "fit_derivative.svg"), lightfig)
save(joinpath(figdir, "example", "fit_derivative.svg"), lightfig)

darkfig = with_theme(_ci, github_dark)
save(joinpath(figdir, "readme", "dark", "fit_derivative.svg"), darkfig)
@info "Saved fit derivative curve"

##
