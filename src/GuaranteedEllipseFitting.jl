module GuaranteedEllipseFitting

using LinearAlgebra
using Parameters
using Random
using StaticArrays
using Setfield

include("util.jl")
include("types.jl")
include("observations.jl")
include("conversion.jl")
include("sampling.jl")
include("normalise_data_context.jl")
include("direct_ellipse_fit.jl")
include("guaranteed_ellipse_fit.jl")
include("fit_ellipse.jl")
include("covariance.jl")

export Observations,
       UncertainObservations,
       NormaliseDataContext,
       ToNormalisedSpace,
       FromNormalisedSpace,
       DirectEllipseFit,
       GuaranteedEllipseFit,
       AlgebraicToGeometric,
       GeometricToAlgebraic,
       LevenbergMarquardt,
       sample_points_on_ellipse,
       fit_ellipse,
       determine_algebraic_covariance,
       determine_geometric_covariance
end
