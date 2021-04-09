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
include("normalize_data_context.jl")
include("direct_ellipse_fit.jl")
include("guaranteed_ellipse_fit.jl")
include("fit_ellipse.jl")
include("covariance.jl")

export Observations,
       NormalizeDataContext,
       ToNormalizedSpace,
       FromNormalizedSpace,
       DirectEllipseFit,
       GuaranteedEllipseFit,
       AlgebraicToGeometric,
       GeometricToAlgebraic,    
       sample_points_on_ellipse,   
       fit_ellipse,
       determine_algebraic_covariance,
       determine_geometric_covariance
end
