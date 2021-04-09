abstract type AbstractContext end
abstract type AbstractNormalizeDataContext <: AbstractContext end
abstract type AbstractNormalizationMethod end
abstract type AbstractCovarianceMatrix end
abstract type AbstractFittingMethod end
abstract type AbstractObservations end
abstract type AbstractUncertainObservations <: AbstractObservations end
abstract type AbstractGauge end
abstract type AbstractOptimizationScheme end
abstract type AbstractLatentEllipseParametrisation end
abstract type AbstractEllipseParameterConversion end

struct AlgebraicToGeometric <: AbstractEllipseParameterConversion end
struct GeometricToAlgebraic <: AbstractEllipseParameterConversion end


struct FirstLatentEllipseParametrisation <: AbstractLatentEllipseParametrisation end
struct SecondLatentEllipseParametrisation <: AbstractLatentEllipseParametrisation end

@with_kw struct ManualEstimation{T <: Union{AbstractVecOrMat, Nothing}} <: AbstractOptimizationScheme
    𝛉::T = nothing
end

struct IsotropicScalingTranslation <: AbstractNormalizationMethod end
struct IsotropicScaling <: AbstractNormalizationMethod end

struct DirectEllipseFit <: AbstractFittingMethod end

@with_kw struct LevenbergMarquardt{T <:  AbstractOptimizationScheme} <: AbstractOptimizationScheme
    max_iter::Integer = 10000
    γ::Float64 = 1.2
    γ_inc::Float64 = 10
    γ_dec::Float64 = 0.1
    λ::Float64 = 0.01
    λ_lower_bound::Float64 = 1e-10
    tol_𝛉::Float64 = 1e-8             # change in parameters tolerance
    tol_cost::Float64 = 1e-12         # change in cost tolerance
    tol_Δ::Float64 = 1e-8             # change in update tolerance
    tol_∇::Float64 = 1e-12            # change in gradient tolerance
    seed::T = ManualEstimation()
end

@with_kw struct GuaranteedEllipseFit <: AbstractFittingMethod 
    parametrisation::AbstractLatentEllipseParametrisation = SecondLatentEllipseParametrisation()
    optimisation_scheme::LevenbergMarquardt = LevenbergMarquardt()
end

GuaranteedEllipseFit(𝛉::AbstractVector) = GuaranteedEllipseFit(optimisation_scheme = LevenbergMarquardt(seed = ManualEstimation(𝛉)))


struct UnitNorm <: AbstractGauge end
struct ToNormalizedSpace end
struct FromNormalizedSpace end



struct JacobianMatrix{T₁ <: Union{Nothing, AbstractObservations}, T₂ <: AbstractMatrix}
    observations::T₁
    jacobian::T₂
end