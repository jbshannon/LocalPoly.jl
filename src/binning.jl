function linear_binning(x, y; nbins=floor(Int, length(x)/100))
    xmin, xmax = extrema(x)
    Δ = abs(xmax - xmin)/(nbins-1) # ensure maximum(x) < maximum(g)
    g = collect(range(xmin, step=Δ, length=nbins))

    Y = zeros(eltype(g), nbins)
    c = zeros(eltype(g), nbins)
    @turbo for i in eachindex(x, y)
        L = (x[i] - xmin)/Δ + 1 # transformation matching gridpoints to indices
        𝓁 = floor(Int, L) # index of left gridpoint
        w = L - 𝓁 # remainder, used for weighting
        Y[𝓁] += w * y[i]
        Y[𝓁+1] += (1-w) * y[i]
        c[𝓁] += w
        c[𝓁+1] += 1-w
    end

    # Replace weighted sum in Y with weighted average
    for i in eachindex(Y, c)
        if c[i] > 0
            Y[i] /= c[i]
        end
    end

    return g, Y, c
end

# TODO: combine different binning methods into one function
function simple_binning(x, y; nbins=floor(Int, length(x)/100))
    xmin, xmax = extrema(x)
    Δ = abs(xmax - xmin)/(nbins - 1) # ensure maximum(x) < maximum(g)
    g = collect(range(xmin, step=Δ, length=nbins))

    Y = zeros(eltype(g), nbins)
    c = zeros(eltype(g), nbins)
    @turbo for i in eachindex(x, y)
        L = (x[i] - xmin)/Δ + 1
        𝓁 = round(Int, L)
        Y[𝓁] += y[i]
        c[𝓁] += 1
    end

    # replace weighted sum in Y with weighted average
    for i in eachindex(Y, c)
        if c[i] > 0
            Y[i] /= c[i]
        end
    end

    return g, Y, c
end
