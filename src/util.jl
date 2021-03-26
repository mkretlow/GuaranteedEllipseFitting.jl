const  √ = sqrt
const  ∑ = sum
const  ⊗ = kron

function hom⁻¹(v::StaticVector)
    if isapprox(v[end], 0.0; atol = 1e-14)
        pop(v)
    else
        pop(v / v[end])
    end
end

function hom(v::StaticVector)
    push(v,1)
end

function svectors(x::AbstractArray{T}, ::Val{N}) where {T,N}
    size(x,1) == N || error("sizes mismatch")
    isbitstype(T) || error("use for bitstypes only")
    reinterpret(SVector{N,T}, vec(x))
end