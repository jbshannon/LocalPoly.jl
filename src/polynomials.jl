"""
    allmultiexponents(N, degree)

All results of `Combinatorics.multiexponents(N, p)` for p ≤ degree.

## Examples
```julia-repl
julia> multiexponents(2, 0) |> collect
1-element Vector{Any}:
 [0, 0]

julia> multiexponents(2, 1) |> collect
2-element Vector{Any}:
 [1, 0]
 [0, 1]

julia> multiexponents(2, 2) |> collect
3-element Vector{Any}:
 [2, 0]
 [1, 1]
 [0, 2]

julia> allmultiexponents(2, 2)
6-element Vector{Tuple{Int64, Int64}}:
 (0, 0)
 (1, 0)
 (0, 1)
 (2, 0)
 (1, 1)
 (0, 2)
````
"""
function allmultiexponents(N, degree)
    multexp = NTuple{N, Int64}[]
    for k in 0:degree
        for mex in multiexponents(N, k)
            push!(multexp, Tuple(mex))
        end
    end
    return multexp
end

function matchexponents(a, b, c)
    return all(i -> a[i] + b[i] == c[i], eachindex(a, b, c))
end

function findpredecessors(E, i)
    for j in 2:i-1, k in 2:j
        matchexponents(E[j], E[k], E[i]) && return (i, j, k)
    end
end

function polypredecessors(E)
    istart = findfirst(>(1) ∘ sum, E)
    inds = Iterators.drop(eachindex(E), istart-1)
    return map(Base.Fix1(findpredecessors, E), inds)
end

function _polybasis!(X, x, x₀)
    for i in eachindex(x)
        x̂ = x[i] - x₀
        x̂ᵖ = 1
        for p in 1:size(X, 2)-1
            x̂ᵖ *= x̂
            X[i, p+1] = x̂ᵖ
        end
    end
    return X
end

function _polybasis(x, x₀, degree)
    X = ones(eltype(x), length(x), degree+1)
    _polybasis!(X, x, x₀)
    return X
end
