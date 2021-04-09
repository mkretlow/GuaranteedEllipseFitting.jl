using GuaranteedEllipseFitting
using Test, StaticArrays, LinearAlgebra

@testset "GuaranteedEllipseFitting.jl" begin
    include("data_points.jl")
    include("test_sampling.jl")
    include("test_observations.jl")
    include("test_conversion.jl")
    include("test_normalize_data_context.jl")
    include("test_direct_ellipse_fit.jl")
    include("guaranteed_ellipse_fit.jl")
    include("test_covariance.jl")
end
