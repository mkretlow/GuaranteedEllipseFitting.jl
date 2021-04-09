@testset "Test Sampling" begin
    @info "Test: Sampling Points on Ellipse"

    @testset "Samples points satisfy ellipse equation" begin
        𝛏 = [10, 5, 25, 25, π/4]
        geo_to_alg = GeometricToAlgebraic()
        𝛉 = geo_to_alg(𝛏)
        α₁ = 0
        α₂ = 2*π
        N = 25
        observations = sample_points_on_ellipse(𝛏..., N, α₁, α₂)
        𝐫 = GuaranteedEllipseFitting.vector_valued_objective(observations, 𝛉)
        @test isapprox(dot(𝐫, 𝐫), 0.0, atol= 1e-10) 
    end 

end