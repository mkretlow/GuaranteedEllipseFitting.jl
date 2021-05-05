abstract type AbstractContext end
abstract type AbstractNormaliseDataContext <: AbstractContext end
abstract type AbstractNormalisationMethod end
abstract type AbstractCovarianceMatrix end
abstract type AbstractFittingMethod end
abstract type AbstractObservations end
abstract type AbstractUncertainObservations <: AbstractObservations end
abstract type AbstractGauge end
abstract type AbstractOptimisationScheme end
abstract type AbstractLatentEllipseParametrisation end
abstract type AbstractEllipseParameterConversion end
abstract type AbstractOptimisationResult end

@with_kw struct OptimisationResult{T₁ <: AbstractVector} <: AbstractOptimisationResult
    minimiser::T₁ = SVector(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    minimum::Float64 = 0.0
    info_str::String = ""
end

struct AlgebraicToGeometric <: AbstractEllipseParameterConversion end
struct GeometricToAlgebraic <: AbstractEllipseParameterConversion end


struct FirstLatentEllipseParametrisation <: AbstractLatentEllipseParametrisation end
struct SecondLatentEllipseParametrisation <: AbstractLatentEllipseParametrisation end

@with_kw struct ManualEstimation{T <: Union{AbstractVecOrMat, Nothing}} <: AbstractOptimisationScheme
    𝛉::T = nothing
end

struct IsotropicScalingTranslation <: AbstractNormalisationMethod end
struct IsotropicScaling <: AbstractNormalisationMethod end

struct DirectEllipseFit <: AbstractFittingMethod end

@with_kw struct LevenbergMarquardt{T <:  AbstractOptimisationScheme} <: AbstractOptimisationScheme
    max_iter::Integer = 10000
    max_func_eval::Integer = 500
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

"""
```
    GuaranteedEllipseFit(; parametrisation = SecondLatentEllipseParametrisation(), optimisation_scheme = LevenbergMarquardt())
```

Configures the choice of implicit parametrisation and the optimisation algorithm; currently only the default values are supported.
You can also choose to display summary information about the optimisation process in the console by setting `verbose = true`.

"""
@with_kw struct GuaranteedEllipseFit <: AbstractFittingMethod
    parametrisation::AbstractLatentEllipseParametrisation = SecondLatentEllipseParametrisation()
    optimisation_scheme::LevenbergMarquardt = LevenbergMarquardt()
    verbose::Bool = false
end

GuaranteedEllipseFit(𝛉::AbstractVector) = GuaranteedEllipseFit(optimisation_scheme = LevenbergMarquardt(seed = ManualEstimation(𝛉)))

struct UnitNorm <: AbstractGauge end
struct ToNormalisedSpace end
struct FromNormalisedSpace end

struct JacobianMatrix{T₁ <: Union{Nothing, AbstractObservations}, T₂ <: AbstractMatrix}
    observations::T₁
    jacobian::T₂
end
