@testset "Test Normalisation" begin
    @info "Test: Ellipse Parameter Cooordinate System Conversion"

    @testset "From Unnormalised to Normalised and back" begin
        𝛏₀ = [10, 5, 25, 25, π/4]
        geo_to_alg = GeometricToAlgebraic()
        s = 0.01
        m₁ = 25
        m₂ = 25
        𝐓 = SMatrix{3,3,Float64}(s,0,0,0,s,0,-s*m₁,-s*m₂,1)
        normalise = NormalizeDataContext(tuple(𝐓))        
        𝛉₀ = geo_to_alg(𝛏₀)
        𝛉₀ = 𝛉₀ / norm(𝛉₀)
        𝛉′ = normalise(ToNormalizedSpace(), 𝛉₀)
        𝛉₁ = normalise(FromNormalizedSpace(), 𝛉′)
        @test all(𝛉₀ .≈ 𝛉₁)
    end
end

