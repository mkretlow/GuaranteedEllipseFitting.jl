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

function transform_covariance(Λ::AbstractVector, 𝐓::AbstractMatrix, ndim::Val{2})
    Λ₂ = map(Λ) do 𝚲
       # Lift the covariance matrix so that it correspond to homeogenous coordinates.
       # This way the requisite transformation can be computed by multiply with a matrix 𝐓.
       𝚲₀ =  hcat(𝚲, SVector(0.0, 0.0))
       𝚲₁ = vcat(𝚲₀, transpose(SVector(0.0, 0.0, 0.0)))
       𝚲₂ =  (𝐓 * 𝚲₁ * 𝐓')
       𝚲′ = SMatrix{2,2,Float64,4}(𝚲₂[1], 𝚲₂[2], 𝚲₂[4], 𝚲₂[5])
    end
end

function (coordinate_transformation::NormalizeDataContext)(direction::FromNormalizedSpace, 𝛉′::AbstractVector)
    𝒯 = matrices(coordinate_transformation)
    𝐓 = 𝒯[1]
    𝐄 = Diagonal(SVector(1, 2^-1, 1, 2^-1, 2^-1, 1))
    # Permutation matrix for interchanging the 3rd and 4th entries of a length-6 vector. 
    𝐏₃₄ = Diagonal(SVector(0,1,0)) ⊗ SMatrix{2,2,Float64}(0,1,1,0) + Diagonal(SVector(1,0,1)) ⊗ SMatrix{2,2,Float64}(1,0,0,1)
    # 9 x 6 duplication matrix
    𝐃₃ = [1 0 0 0 0 0; 
          0 1 0 0 0 0;
          0 0 1 0 0 0;
          0 1 0 0 0 0;
          0 0 0 1 0 0;
          0 0 0 0 1 0;
          0 0 1 0 0 0;
          0 0 0 0 1 0;
          0 0 0 0 0 1]    
    𝛉 = 𝐄 \ 𝐏₃₄ * pinv(𝐃₃) * kron(𝐓, 𝐓)' * 𝐃₃ * 𝐏₃₄ * 𝐄 * 𝛉′
    𝛉 = 𝛉 / norm(𝛉)
    return 𝛉   
end

function (coordinate_transformation::NormalizeDataContext)(direction::ToNormalizedSpace, 𝛉::AbstractVector)
    𝒯 = matrices(coordinate_transformation)
    𝐓 = 𝒯[1]
    𝐄 = Diagonal(SVector(1, 2^-1, 1, 2^-1, 2^-1, 1))
    # Permutation matrix for interchanging the 3rd and 4th entries of a length-6 vector. 
    𝐏₃₄ = Diagonal(SVector(0,1,0)) ⊗ SMatrix{2,2,Float64}(0,1,1,0) + Diagonal(SVector(1,0,1)) ⊗ SMatrix{2,2,Float64}(1,0,0,1)
    # 9 x 6 duplication matrix
    𝐃₃ = [1 0 0 0 0 0; 
          0 1 0 0 0 0;
          0 0 1 0 0 0;
          0 1 0 0 0 0;
          0 0 0 1 0 0;
          0 0 0 0 1 0;
          0 0 1 0 0 0;
          0 0 0 0 1 0;
          0 0 0 0 0 1]    
    𝛉′ = 𝐄 \ 𝐏₃₄ * pinv(𝐃₃) * inv(kron(𝐓, 𝐓))' * 𝐃₃ * 𝐏₃₄ * 𝐄 * 𝛉
    𝛉′ = 𝛉′ / norm(𝛉′)
    return 𝛉′   
end