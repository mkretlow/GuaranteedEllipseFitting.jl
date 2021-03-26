
abstract type AbstractContext end
abstract type AbstractNormalizeDataContext <: AbstractContext end
abstract type AbstractNormalizationMethod end
abstract type AbstractCovarianceMatrix end
abstract type AbstractFittingMethod end
abstract type AbstractObservations end
abstract type AbstractUncertainObservations <: AbstractObservations end
abstract type AbstractGauge end

struct IsotropicScalingTranslation <: AbstractNormalizationMethod end
struct IsotropicScaling <: AbstractNormalizationMethod end

struct DirectEllipseFit <: AbstractFittingMethod end
struct GuaranteedEllipseFit <: AbstractFittingMethod end

struct UnitNorm <: AbstractGauge end
struct ToNormalizedSpace end
struct FromNormalizedSpace end