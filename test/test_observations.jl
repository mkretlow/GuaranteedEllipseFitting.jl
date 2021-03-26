@testset "Observations" begin
    @info "Test: Observations"
    
    @testset "API" begin
        @inferred Observations(tuple(data))
    end
end
