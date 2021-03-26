using GuaranteedEllipseFitting
using Test, StaticArrays

@testset "GuaranteedEllipseFitting.jl" begin
    include("data_points.jl")
    include("test_observations.jl")
    #include("test_normalize_data_context.jl")
    include("test_direct_ellipse_fit.jl")
end
