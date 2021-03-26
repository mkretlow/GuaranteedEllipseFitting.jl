module GuaranteedEllipseFitting

using LinearAlgebra
using Parameters
using Random
using StaticArrays

include("util.jl")
include("types.jl")
include("observations.jl")
include("normalize_data_context.jl")
include("direct_ellipse_fit.jl")
include("fit_ellipse.jl")

export Observations,
       NormalizeDataContext,
       DirectEllipseFit,
       GuaranteedEllipseFit,
       fit_ellipse
end
