@testset "Covariance Matrix for GuaranteedEllipseFit" begin
    @info "Test: Guaranteed Ellipse Fit Covariance Matrix"

    @testset "Numerical" begin
        𝚺_ref = [9.26328248166718          -2.0798069234395        -0.583656167803226        0.0481226621171429        0.0018586684115733;
                -2.07980692343951            6.006823008749        0.0378537156164336       0.00295257680184317     -0.000656827363252417;
                -0.583656167803228        0.0378537156164377          7.37460418632127         0.158598345149172      0.000644886416916763;
                 0.0481226621171417       0.00295257680184355         0.158598345149172          3.72310516565684       -0.0007869916229663;
                 0.00185866841157333     -0.000656827363252431      0.000644886416916745     -0.000786991622966293      0.000179354045153511]


        ℳ = GuaranteedEllipseFitting.svectors(transpose(data), Val{2}())
        observations = Observations(tuple(ℳ))

        𝛉₀ = fit_ellipse(observations, DirectEllipseFit())
        𝛉₁ = fit_ellipse(observations, GuaranteedEllipseFit(𝛉₀))
        alg_to_geo = AlgebraicToGeometric()
        𝛈 = alg_to_geo(𝛉₁)
        𝚲 = determine_algebraic_covariance(observations, 𝛉₁, GuaranteedEllipseFit())
        𝚺 = determine_geometric_covariance(observations, 𝛈, GuaranteedEllipseFit())
        @test isapprox(norm(𝚺_ref - 𝚺), 0.0, atol= 1e-5) 
    end

end
