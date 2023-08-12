"""
    polybasis_grid(L, δ, p)

Construct a polynomial basis of `L` steps of length `delta` and degree `p`.

``A[i, j] = ((j - (size(A, 2) ÷ 2 + 1)) * δ)^(i-1)``

## Examples
```julia-repl
julia> polybasis_grid(3, 0.5, 3)
4×7 Matrix{Float64}:
  1.0     1.0   1.0     1.0  1.0    1.0  1.0
 -1.5    -1.0  -0.5    -0.0  0.5    1.0  1.5
  2.25    1.0   0.25    0.0  0.25   1.0  2.25
 -3.375  -1.0  -0.125  -0.0  0.125  1.0  3.375
```
"""
function polybasis_grid(L, δ, p)
    A = ones(typeof(δ), 1+p, 1+2L)
    for l in 0:L
        dist = 1
        step = l*δ
        for i in 1:p
            dist *= step
            A[1+i, 1+L+l] = dist
            A[1+i, 1+L-l] = (-1)^i * dist
        end
    end
    return A
end

"""
    gridsteps(g, h, τ)

Compute the number of grid steps for which kernel weights are non-zero. `g` should be a `StepRange`.

``L = min { ⌊τδ/h⌋, M-1 }``
"""
function gridsteps(g, h, τ)
    M = length(g)
    δ = step(g)
    L = min(floor(Int, τ*h/δ), M-1)
    return L
end

"""
    outerproduct!(A::Array{T, N}, xs) where {T, N}

Compute the `N`-dimensional outer product of the vectors `xs`, modifying A in-place.

## Examples

```julia-repl
julia> xs = ([0, 1, 2], [1, 2, 3])
([0, 1, 2], [1, 2, 3])

julia> A = xs[1] * xs[2]'
3×3 Matrix{Int64}:
 0  0  0
 1  2  3
 2  4  6

julia> fill!(A, 0)
3×3 Matrix{Int64}:
 0  0  0
 0  0  0
 0  0  0

julia> outerproduct!(A, xs)
3×3 Matrix{Int64}:
 0  0  0
 1  2  3
 2  4  6
```
"""
@generated function outerproduct!(A::Array{T, N}, xs) where {T, N}
    quote
        @nloops $N i A begin
            s = one(T)
            @nexprs $N j -> s *= xs[j][i_j]
            (@nref $N A i) = s
        end
        return A
    end
end

"""
    outerproduct(xs::NTuple{N, T}) where {N, T}

Compute the `N`-dimensional outer product of the vectors `xs`

## Examples

```julia-repl
julia> xs = ([0, 1, 2], [1, 2, 3])
([0, 1, 2], [1, 2, 3])

julia> xs[1] * xs[2]'
3×3 Matrix{Int64}:
 0  0  0
 1  2  3
 2  4  6

julia> outerproduct(xs)
3×3 Matrix{Int64}:
 0  0  0
 1  2  3
 2  4  6
```
"""
function outerproduct(x::NTuple{N, S}) where {N, S}
    A = zeros(eltype(S), length.(x))
    return outerproduct!(A, x)
end

function fillκ!(κ::Array{T, N}, A, kernelweights::Array{T, N}, k) where {T, N}
    Arows = ntuple(i -> view(A[i], 1+k[i], :), Val(N))
    outerproduct!(κ, Arows) # construct grid of terms (lδ)ᵏ
    broadcast!(*, κ, κ, kernelweights) # multiply by kernel weights Kₕ(lδ)
    return κ
end

function fillκ(A, kernelweights, k)
    κ = zero(kernelweights)
    fillκ!(κ, A, kernelweights, k)
end

## Direct convolution

struct ConvolutionData{T <: Real, N}
    grid::GridData
    A::NTuple{N, Matrix{T}}
    degree::Int
    multexp::Vector{NTuple{N, Int}}
    kernels::Vector{Array{T, N}}
    s̃
    t̃
end

function ConvolutionData(
    grid::GridData{T, N, R},
    h;
    degree=1,
    kernel=:Epanechnikov,
) where {T, N, R}
    @unpack g, c, d = grid
    M = length.(g)
    τ = support(Val(kernel))
    L = gridsteps.(g, h, τ)
    A = ntuple(i -> polybasis_grid(L[i], step(g[i]), max(2, 2degree)), Val(N))
    Arows = ntuple(i -> view(A[i], 2, :), N)
    kernelweights = map(u -> Kₕ(Val(kernel), u, h), gridnodes(Arows))

    multexp = allmultiexponents(N, 2degree)
    kernels = map(k -> fillκ(A, kernelweights, k), multexp)
    s̃ = zeros(T, M..., length(multexp))
    t̃ = zeros(T, M..., count(<=(degree) ∘ sum, multexp))
    return ConvolutionData(grid, A, degree, multexp, kernels, s̃, t̃)
end

function show(io::IO, D::ConvolutionData{T, N}) where {T, N}
    @unpack kernels, multexp = D
    degree = sum(last(multexp))÷2
    L = gridsteps(kernels)
    println(io, typeof(D))
    println(io, "      Degree: $degree")
    println(io, "  Grid steps: $(length(L) > 1 ? L : first(L))")
end

gridsteps(kernel) = ntuple(i -> (size(kernel, i)-1)÷2, ndims(kernel))
gridsteps(kernels::NTuple{N, T}) where {N, T} = gridsteps(first(kernels))
gridsteps(C::ConvolutionData) = gridsteps(first(C.kernels))

function convolve!(C, A, B)
    Bmid = CartesianIndices(B)[length(B)÷2+1]
    for I in CartesianIndices(A), J in CartesianIndices(B)
        I′ = I + J - Bmid
        if I′ in CartesianIndices(A)
            C[I] += A[I′] * B[J]
        end
    end
end

function convolve!(convdata::ConvolutionData{T, N}) where {T, N}
    @unpack grid, degree, multexp, kernels, s̃, t̃ = convdata
    @unpack c, d = grid
    colons = ntuple(i -> Base.Colon(), N)
    Threads.@threads for i in eachindex(kernels)
        convolve!(view(s̃, colons..., i), c, kernels[i])
        if sum(multexp[i]) <= degree
            convolve!(view(t̃, colons..., i), d, kernels[i])
        end
    end
    return convdata
end

function finalstep(s̃, t̃, multexp)
    p = ndims(t̃) == 1 ? 1 : last(size(t̃))
    MM = [multexp[i] .+ multexp[j] for j in 1:p, i in 1:p]
    MI = map(m -> findfirst(==(m), multexp), MM)
    XWX = MMatrix{p, p}(zeros(eltype(s̃), p, p))
    colons = ntuple(i -> Base.Colon(), length(first(multexp)))
    map(CartesianIndices(view(s̃, colons..., 1))) do I
        copyto!(XWX, view(s̃, Tuple(I)..., MI))
        return lu(XWX) \ view(t̃, Tuple(I)..., :)
    end
end

function finalstep(data::ConvolutionData)
    @unpack s̃, t̃, multexp = data
    finalstep(s̃, t̃, multexp)
end
