 @testset "Direct Ellipse Fit" begin
    @info "Test: Direct Ellipse Fit"

    @testset "Numerical" begin
        truth = SVector(4.053658735425687e-6,
                        -1.3487952962436518e-7,
                        1.8534714216409155e-5,
                        -0.001965777847001004,
                        -0.009242570343232327,
                         0.9999553541288334)

        𝛉 = fit_ellipse(data, DirectEllipseFit())
        @test all(𝛉 .≈ truth)
    end

end

