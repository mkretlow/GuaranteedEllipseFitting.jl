struct NormalizeDataContext{N,T <:AbstractMatrix} <: AbstractNormalizeDataContext
    transformations::NTuple{N,T}
end

function matrices(context::NormalizeDataContext)
    context.transformations
end

NormalizeDataContext(observations::AbstractObservations, method::AbstractNormalizationMethod) = NormalizeDataContext(method(observations))

function (method::AbstractNormalizationMethod)(observations::AbstractObservations)
    return method(observations)
end

function (method::IsotropicScalingTranslation)(observations::AbstractObservations)
    @unpack data = observations
    𝐓 = isotropic_scale_and_translate(data[1])
    return tuple(𝐓)
end

function (method::IsotropicScaling)(observations::AbstractObservations)
    @unpack data = observations
    𝐓 = isotropic_scale(data[1])
    return tuple(𝐓)
end

function (normalize_data::NormalizeDataContext)(observations::Observations) 
    transformations = matrices(normalize_data)
    data = _normalize_data(observations, transformations)
    return Observations(data)
end


function (normalize_data::NormalizeDataContext)(observations::UncertainObservations) 
    transformations = matrices(normalize_data)
    data = _normalize_data(observations, transformations)
    covariance_matrices = _normalize_covariance_matrices(observations, transformations)
    return UncertainObservations(data, covariance_matrices)
end

function _normalize_data(observations::AbstractObservations, transformations)
    @unpack data = observations
    𝐓₁ = transformations[1]
    ℳ = data[1]
    return tuple(transform_data(ℳ, 𝐓₁))
end

function  _normalize_covariance_matrices(observations::AbstractObservations, transformations)
    @unpack covariance_matrices = observations
    𝐓₁ = transformations[1]
    Λ = covariance_matrices[1]
    D = size(first(Λ),1)
    return tuple(transform_covariance(Λ, 𝐓₁, Val(D)))
end


function isotropic_scale_and_translate(ℳ::AbstractVector)
    if isempty(ℳ)
        throw(ArgumentError("Array cannot be empty."))
    end
    npts = length(ℳ)
    ndim = length(ℳ[1])
    𝐜 = centroid(ℳ)
    σ = root_mean_square(ℳ, 𝐜)
    σ⁻¹ = 1 / σ
    𝐓 = SMatrix{ndim+1,ndim+1,Float64, (ndim+1)^2}([σ⁻¹*Matrix{Float64}(I,ndim,ndim) -σ⁻¹*𝐜 ; zeros(1,ndim) 1.0])
end

function isotropic_scale(ℳ::AbstractVector)
    if isempty(ℳ)
        throw(ArgumentError("Array cannot be empty."))
    end
    npts = length(ℳ)
    ndim = length(ℳ[1])
    𝐜 = centroid(ℳ)
    σ = root_mean_square(ℳ, 𝐜)
    σ⁻¹ = 1 / σ
    𝐓 = SMatrix{ndim+1,ndim+1,Float64, (ndim+1)^2}([σ⁻¹*Matrix{Float64}(I,ndim,ndim) zeros(ndim) ; zeros(1,ndim) 1.0])
end

function centroid(positions::AbstractVector) 
    x = zeros(eltype(positions))
    for pos ∈ positions
        x = x + pos
    end
    return x / length(positions)
end

function root_mean_square(ℳ::AbstractVector, 𝐜::AbstractVector) 
    total = 0.0
    npts = length(ℳ)
    ndim = length(ℳ[1])
    for 𝐦 ∈ ℳ
         total  = total + ∑((𝐦-𝐜).^2)
    end
    σ = √( (1/(npts*ndim)) * total)
end

function transform_data(ℳ::AbstractVector, 𝐓::AbstractMatrix)
    𝒪 = transform_data!(copy(ℳ), 𝐓)
end

function transform_data!(ℳ::AbstractVector, 𝐓::AbstractMatrix)
    map!(ℳ , ℳ) do 𝐦
         hom⁻¹(𝐓 * hom(𝐦))
    end
     ℳ
end

function transform_data(ℳ::AbstractVector)
    𝒪, 𝐓 = transform_data!(copy(ℳ))
end

function transform_data!(ℳ::AbstractVector)
    𝐓 = transform_data(ℳ)
    map!(ℳ , ℳ) do 𝐦
         hom⁻¹(𝐓 * hom(𝐦))
    end
     ℳ, 𝐓
end


function (coordinate_transformation::NormalizeDataContext)(direction::FromNormalizedSpace, 𝛉::AbstractVector)
    𝒯 = matrices(coordinate_transformation)
    𝐓 = 𝒯[1]

    # TODO replace this with a more elegant (direct) formula.
    a, b, c, d, e, f = 𝛉
    𝐂 =[a b/2 d/2;  b/2 c e/2 ; d/2 e/2 f]
    𝐂′ = 𝐓' * 𝐂 * 𝐓
    𝛉′ = SVector(𝐂′[1,1], 𝐂′[1,2]*2, 𝐂′[2,2], 𝐂′[1,3]*2, 𝐂′[2,3]*2, 𝐂′[3,3])
    return 𝛉′ / norm(𝛉′)   
end