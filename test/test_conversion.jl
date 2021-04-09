@testset "Test Conversion" begin
    @info "Test: Ellipse Parameter Conversion"

    @testset "From Geometric to Algebraic and back" begin
        𝛏₀ = [10, 5, 25, 25, π/4]
        geo_to_alg = GeometricToAlgebraic()
        alg_to_geo = AlgebraicToGeometric()
        𝛉 = geo_to_alg(𝛏₀)
        𝛏₁ = alg_to_geo(𝛉)
        @test all(𝛏₀ .≈ 𝛏₁)
    end

    @testset "From Algebraic to Geometric and back" begin
        𝛉₀ = [0.025, -0.03, 0.025, -0.5, -0.5, 11.5]
        geo_to_alg = GeometricToAlgebraic()
        alg_to_geo = AlgebraicToGeometric()
        𝛏 = alg_to_geo(𝛉₀)
        𝛉₁ = geo_to_alg(𝛏)
        @test all(𝛉₀ .≈ 𝛉₁)
    end

    @testset "From  Algebraic to Latent Algebraic and back" begin
        𝛏₀ = [10, 5, 25, 25, π/4]
        geo_to_alg = GeometricToAlgebraic()
        𝛉₀ = geo_to_alg(𝛏₀)
        𝛉₀ = 𝛉₀ / norm(𝛉₀)

        𝛋 = GuaranteedEllipseFitting.to_latent_parameters(GuaranteedEllipseFitting.FirstLatentEllipseParametrisation(), 𝛉₀)
        𝛉₁ = GuaranteedEllipseFitting.from_latent_parameters(GuaranteedEllipseFitting.FirstLatentEllipseParametrisation(), 𝛋)
        @test all(𝛉₀ .≈ 𝛉₁)

        𝛋 = GuaranteedEllipseFitting.to_latent_parameters(GuaranteedEllipseFitting.SecondLatentEllipseParametrisation(), 𝛉₀)
        𝛉₁ = GuaranteedEllipseFitting.from_latent_parameters(GuaranteedEllipseFitting.SecondLatentEllipseParametrisation(), 𝛋)
        @test all(𝛉₀ .≈ 𝛉₁)
    end

end