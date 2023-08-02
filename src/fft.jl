

# copied from https://oeis.org/A002182/list
const HIGHLY_COMPOSITE_NUMBERS = [
    1,2,4,6,12,24,36,48,60,120,180,240,360,720,840,
    1260,1680,2520,5040,7560,10080,15120,20160,25200,
    27720,45360,50400,55440,83160,110880,166320,
    221760,277200,332640,498960,554400,665280,720720,
    1081080,1441440,2162160
]

## FFT convolutions

struct FFTData{T, N, C <: Complex}
    data::ConvolutionData{T, N}
    cᶻ::Array{T, N}
    dᶻ::Array{T, N}
    tmp::Array{T, N}
    ℱ # FFT plan
    Cᶻ::Array{C, N}
    Dᶻ::Array{C, N}
    Kᶻ
    TMP::Array{C, N}
end

getfirst(f, x) = x[findfirst(f, x)]

"""
    padκ!(κᶻ::Array{T, N}, κ::Array{T, N}, tmp::Array{T, N}) where {T, N}

Zero-pad and shift the convolution kernel, modifying `κᶻ` in-place.
"""
function padκ!(κᶻ::Array{T, N}, κ::Array{T, N}, tmp::Array{T, N}) where {T, N}
    copyto!(view(tmp, CartesianIndices(κ)), reverse!(κ)) # copy reversed kernel matrix
    circshift!(κᶻ, tmp, ntuple(i -> -(size(κ, i) - 1)/2, Val(N))) # shift by -L
    return κᶻ
end

"""
    padκ(κ::Array{T, N}, tmp::Array{T, N}) where {T, N}

Zero-pad and shift the convolution kernel.
"""
function padκ(κ::Array{T, N}, tmp::Array{T, N}) where {T, N}
    κᶻ = zero(tmp)
    return padκ!(κᶻ, κ, tmp)
end

function FFTData(data::ConvolutionData{T, N}) where {T, N}
    @unpack grid, kernels = data
    @unpack g, c, d = grid

    M = length.(g)
    L = gridsteps(kernels)
    P = ntuple(i -> getfirst(>=(M[i] + L[i]), HIGHLY_COMPOSITE_NUMBERS), Val(N))
    cᶻ = zeros(T, P)
    dᶻ = zeros(T, P)
    copyto!(view(cᶻ, CartesianIndices(c)), c)
    copyto!(view(dᶻ, CartesianIndices(d)), d)

    # Pre-allocate FFT arrays
    ℱ = plan_rfft(cᶻ)
    Cᶻ = ℱ * cᶻ
    Dᶻ = ℱ * dᶻ

    # Pre-compute FFT of kernels
    tmp = zero(cᶻ)
    TMP = zero(Cᶻ)
    Kᶻ = map(κ -> ℱ * padκ(κ, tmp), kernels)
    return FFTData(data, cᶻ, dᶻ, tmp, ℱ, Cᶻ, Dᶻ, Kᶻ, TMP)
end

function show(io::IO, D::FFTData{T, N, C}) where {T, N, C}
    @unpack Cᶻ = D
    P = ndims(Cᶻ) > 1 ? size(Cᶻ) : length(Cᶻ)
    println(io, typeof(D))
    println(io, " Size: $P")
end

function FFTData(
    y, x, h;
    nbins=guessbins(x),
    degree=1,
    kernel=:Epanechnikov,
)
    grid = linear_binning(y, x; nbins)
    return FFTData(ConvolutionData(grid, h; degree, kernel))
end

function convolve!(dest, A, B, C, plan, c)
    broadcast!(*, C, A, B) # point-wise multiply the two matrices
    ldiv!(c, plan, C) # inverse FFT
    copyto!(dest, view(c, CartesianIndices(dest))) # extract non-zero-padded entries
end

function convolve!(fftdata::FFTData)
    @unpack data, Kᶻ, Cᶻ, Dᶻ, ℱ, TMP, tmp = fftdata
    @unpack degree, multexp, s̃, t̃ = data
    colons = ntuple(i -> Base.Colon(), ndims(s̃)-1)
    for i in eachindex(Kᶻ) # threading doesn't add much improvement
        convolve!(view(s̃, colons..., i), Cᶻ, Kᶻ[i], TMP, ℱ, tmp)
        if sum(multexp[i]) <= degree
            convolve!(view(t̃, colons..., i), Dᶻ, Kᶻ[i], TMP, ℱ, tmp)
        end
    end
    return fftdata
end

function finalstep(FFT::FFTData)
    @unpack data = FFT
    finalstep(data)
end
