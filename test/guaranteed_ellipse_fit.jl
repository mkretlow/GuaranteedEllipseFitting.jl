@testset "Guaranteed Ellipse Fit" begin
    @info "Test: Guaranteed Ellipse Fit"

    @testset "Convergence on noiseless data" begin
        𝛏 = [10, 5, 25, 25, π/4]
        geo_to_alg = GeometricToAlgebraic()
        𝛉₀ = geo_to_alg(𝛏)
        𝛉₀ = 𝛉₀ / norm(𝛉₀)
        α₁ = 0
        α₂ = 2*π
        N = 25
        observations = sample_points_on_ellipse(𝛏..., N, α₁, α₂)
        𝛉₁ = fit_ellipse(observations, GuaranteedEllipseFit(𝛉₀))
        𝐫 = GuaranteedEllipseFitting.vector_valued_objective(observations, 𝛉₁)
        @test 𝛉₀ ≈ 𝛉₁        
        @test isapprox(dot(𝐫, 𝐫), 0.0, atol= 1e-10) 
    end

    @testset "Numerical" begin
        truth = SVector(3.96548501597781e-06,
                        -5.0179230231193e-07,
                        1.88176818319181e-05,
                        -0.00182736162233016,
                        -0.00930825073431271,
                         0.999955007411677)
        truth_geo = SVector(314.346509735058,
                            144.209111864029,
                            246.264470570879,
                            250.610687177739,
                            0.0168864414911359)

        𝛉₀ = fit_ellipse(data, DirectEllipseFit())
        𝛉₁ = fit_ellipse(data, GuaranteedEllipseFit(𝛉₀))
        alg_to_geo = AlgebraicToGeometric()
        𝛈 = alg_to_geo(𝛉₁)
        @test isapprox(norm(𝛉₁ - truth) , 0.0, atol= 1e-8) 
        @test isapprox(norm(𝛈 - truth_geo) , 0.0, atol= 1e-4) 
    end


end
