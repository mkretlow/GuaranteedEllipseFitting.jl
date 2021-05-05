function estimate_noise_level(observations::Observations, 𝛉::AbstractVector)
    @unpack data = observations
    𝛉 = 𝛉 / norm(𝛉)
    ℳ = data[1]
    N = length(ℳ)
    aml = 0.0
    for n = 1:N
        𝐦 = ℳ[n]
        𝐮ₙ = SA_F64[𝐦[1]^2, 𝐦[1]*𝐦[2], 𝐦[2]^2, 𝐦[1], 𝐦[2], 1]
        ∂𝐮ₙ = SA_F64[2*𝐦[1]  𝐦[2]  0  1  0  0; 0 𝐦[1] 2*𝐦[2] 0 1 0]'
        𝐀ₙ = 𝐮ₙ * 𝐮ₙ'
        𝚲ₙ =  SA_F64[1 0 ; 0 1]
        𝐁ₙ = ∂𝐮ₙ * 𝚲ₙ * ∂𝐮ₙ'
        aml = aml + abs((𝛉' * 𝐀ₙ * 𝛉)/(𝛉' * 𝐁ₙ * 𝛉))
    end
    σ² = aml / (N-5)
    return σ²
end

function determine_algebraic_covariance(observations::Observations, 𝛉₀::AbstractVector, estimator::GuaranteedEllipseFit)
    @unpack data = observations
    # We derive the covariance matrix of the estimate under the unit-norm Gauge constraint.
    𝛉₀ = 𝛉₀ / norm(𝛉₀)
    N = length(first(data))
    σ² = estimate_noise_level(observations, 𝛉₀)
    Λ = [SA_F64[σ² 0.0; 0.0 σ²] for n = 1:N]
    uncertain_observations = UncertainObservations(data, tuple(Λ))
    return determine_algebraic_covariance(uncertain_observations, 𝛉₀, estimator)
end

function determine_algebraic_covariance(observations::UncertainObservations, 𝛉₀::AbstractVector, estimator::GuaranteedEllipseFit)
     # Convert observations and covariance matrices to a data-driven normalised coordinate system.
    normalise = NormaliseDataContext(observations, IsotropicScalingTranslation())
    𝒯 = matrices(normalise)
    𝐓 = 𝒯[1]
    𝛉₁ = normalise(ToNormalisedSpace(), 𝛉₀)
    @unpack data, covariance_matrices = normalise(observations)
    ℳ = data[1]
    Λ = first(covariance_matrices)
    𝐌 = zeros(6,6)
    N = length(first(data))
    for n = 1:N
        𝐦 = ℳ[n]
        𝚲ₙ = Λ[n]
        𝐮ₙ = SA_F64[𝐦[1]^2, 𝐦[1]*𝐦[2], 𝐦[2]^2, 𝐦[1], 𝐦[2], 1]
        ∂𝐮ₙ = SA_F64[2*𝐦[1]  𝐦[2]  0  1  0  0; 0 𝐦[1] 2*𝐦[2] 0 1 0]'
        𝐀ₙ = 𝐮ₙ * 𝐮ₙ'
        𝐁ₙ = ∂𝐮ₙ * 𝚲ₙ * ∂𝐮ₙ'
        𝐌 = 𝐌 + 𝐀ₙ / (𝛉₁' * 𝐁ₙ * 𝛉₁)
    end
    𝐏ₜ = UniformScaling(1) - (𝛉₁*𝛉₁') / norm(𝛉₁)^2
    # Compute rank-5 constrained pseudo-inverse of 𝐌.
    F = svd(𝐌)
    𝐒 = vcat(SVector([1 / F.S[i] for i = 1:5]...), 0)
    𝐌⁻¹ = F.U * Diagonal(𝐒) * F.V'
    𝚺₁ = 𝐏ₜ * 𝐌⁻¹ * 𝐏ₜ

    # Matrices used to transform from normalised to unnormalised (original) coordinate system.
    𝐄 = Diagonal(SVector(1, 2^-1, 1, 2^-1, 2^-1, 1))
    # Permutation matrix for interchanging the 3rd and 4th entries of a length-6 vector.
    𝐏₃₄ = Diagonal(SVector(0,1,0)) ⊗ SMatrix{2,2,Float64}(0,1,1,0) + Diagonal(SVector(1,0,1)) ⊗ SMatrix{2,2,Float64}(1,0,0,1)
    # 9 x 6 duplication matrix
    𝐃₃ = [1 0 0 0 0 0;
          0 1 0 0 0 0;
          0 0 1 0 0 0;
          0 1 0 0 0 0;
          0 0 0 1 0 0;
          0 0 0 0 1 0;
          0 0 1 0 0 0;
          0 0 0 0 1 0;
          0 0 0 0 0 1]
    𝐅 = 𝐄 \ 𝐏₃₄ * pinv(𝐃₃) * kron(𝐓, 𝐓)' * 𝐃₃ * 𝐏₃₄ * 𝐄
    # 𝛉₀ is not unit-normalised
    𝛉₀ = 𝐅 * 𝛉₁
    𝐏ₜ = UniformScaling(1) - (𝛉₀*𝛉₀') / norm(𝛉₀)^2
    𝚺₀ = norm(𝛉₀)^-2 * 𝐏ₜ * 𝐅 * 𝚺₁ * 𝐅' * 𝐏ₜ
    return 𝚺₀
end

function determine_geometric_covariance(observations::Observations, 𝛏::AbstractVector, estimator::GuaranteedEllipseFit)
    geo_to_alg = GeometricToAlgebraic()
    𝛉 =  geo_to_alg(𝛏)
    𝛉 = 𝛉 / norm(𝛉)
    a, b, c, d, e ,f = 𝛉

    Δ = b^2 - 4*a*c
    λ₊ = 0.5*(a + c - (b^2 + (a - c)^2)^0.5)
    λ₋ = 0.5*(a + c + (b^2 + (a - c)^2)^0.5)

    ψ = b*d*e - a*e^2 - b^2*f + c*(4*a*f - d^2)
    V₊ = (ψ/(λ₊*Δ))^0.5
    V₋ = (ψ/(λ₋*Δ))^0.5

    A = max(V₊, V₋)
    B = min(V₊, V₋)

    ∂A = A == V₊ ? ∂V₊(𝛉, ψ, λ₊, Δ) : ∂V₋(𝛉, ψ, λ₋, Δ)
    ∂B = A == V₊ ? ∂V₋(𝛉, ψ, λ₋, Δ) : ∂V₊(𝛉, ψ, λ₊, Δ)

    ∂𝛏 = vcat(∂A', ∂B', ∂H(𝛉, Δ)', ∂K(𝛉, Δ)',  ∂τ(𝛉)')

    𝚺₀ = determine_algebraic_covariance(observations, 𝛉, estimator)
    return ∂𝛏*𝚺₀*∂𝛏'
end

function ∂H(𝛉::AbstractVector, Δ::Number)
    a, b, c, d, e, f = 𝛉
    ∂H_a = (4*c*(2*c*d-b*e))/Δ^2
    ∂H_b = (b^2*e+4*a*c*e-4*b*c*d)/Δ^2
    ∂H_c = (2*b*(b*d-2*a*e))/Δ^2
    ∂H_d = (2*c)/Δ
    ∂H_e = -b/Δ
    ∂H_f =  0
    return SVector(∂H_a, ∂H_b, ∂H_c, ∂H_d, ∂H_e, ∂H_f)
end

function ∂K(𝛉::AbstractVector, Δ::Number)
    a, b, c, d, e, f = 𝛉
    ∂K_a = (2*b*(b*e-2*c*d))/Δ^2
    ∂K_b = (b^2*d+4*a*c*d-4*a*b*e)/Δ^2
    ∂K_c = (4*a*(2*a*e-b*d))/Δ^2;
    ∂K_d = -b/Δ
    ∂K_e =  2*a/Δ
    ∂K_f =  0
    return SVector(∂K_a, ∂K_b, ∂K_c, ∂K_d, ∂K_e, ∂K_f)
end

function ∂τ(𝛉::AbstractVector)
    a, b, c, d, e, f = 𝛉
    ∂τ_a = -b/(2*(b^2+(a-c)^2))
    ∂τ_b = (a-c)/(2*(b^2+(a-c)^2))
    ∂τ_c = b/(2*(b^2+(a-c)^2))
    ∂τ_d = 0
    ∂τ_e = 0
    ∂τ_f = 0
    return SVector(∂τ_a, ∂τ_b, ∂τ_c, ∂τ_d, ∂τ_e, ∂τ_f)
end

function ∂V₊(𝛉::AbstractVector, ψ::Number, λ₊::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    return SVector(∂V₊_a(𝛉, ψ, λ₊, Δ),
                   ∂V₊_b(𝛉, ψ, λ₊, Δ),
                   ∂V₊_c(𝛉, ψ, λ₊, Δ),
                   ∂V₊_d(𝛉, ψ, λ₊, Δ),
                   ∂V₊_e(𝛉, ψ, λ₊, Δ),
                   ∂V₊_f(𝛉, ψ, λ₊, Δ))
end

function ∂V₊_a(𝛉::AbstractVector, ψ::Number, λ₊::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    part1 = 1/(2*λ₊*Δ);
    part2 = (ψ/(λ₊*Δ))^(-0.5);
    part3 = 4*c*f - e^2 + 4*Δ^(-1)*c*ψ;
    part4 = (ψ/(2*λ₊))*(1 + ((c - a)/((a - c)^2 + b^2)^0.5));
    return part1*part2*(part3 - part4)
end

function ∂V₊_b(𝛉::AbstractVector, ψ::Number, λ₊::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    part1 = 1/(2*λ₊*Δ)
    part2 = (ψ/(λ₊*Δ))^(-0.5)
    part3 = d*e - 2*b*f - 2*Δ^(-1)*b*ψ
    part4 = (b*ψ)/((2*λ₊)*((a - c)^2 + b^2)^0.5)
    return  part1*part2*(part3 + part4)
end

function ∂V₊_c(𝛉::AbstractVector, ψ::Number, λ₊::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    part1 = 1/(2*λ₊*Δ)
    part2 = (ψ/(λ₊*Δ))^(-0.5)
    part3 = 4*a*f - d^2 + 4*Δ^(-1)*a*ψ
    part4 = (ψ/(2*λ₊))*(1 + ((a - c)/((a - c)^2 + b^2)^0.5))
    return part1*part2*(part3 - part4)
end

function ∂V₊_d(𝛉::AbstractVector, ψ::Number, λ₊::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    part1 = (ψ/(λ₊*Δ))^(0.5)
    part2 = (b*e - 2*c*d)/(2*ψ)
    return part1*part2;
end

function ∂V₊_e(𝛉::AbstractVector, ψ::Number, λ₊::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    part1 = (ψ/(λ₊*Δ))^(0.5)
    part2 = (b*d - 2*a*e)/(2*ψ)
    return part1*part2;
end

function ∂V₊_f(𝛉::AbstractVector, ψ::Number, λ₊::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    part1 = -(2*λ₊)^-1
    part2 = (ψ/(λ₊*Δ))^(-0.5)
    return part1*part2
end

function ∂V₋(𝛉::AbstractVector, ψ::Number, λ₋::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    return SVector(∂V₋_a(𝛉, ψ, λ₋, Δ),
                   ∂V₋_b(𝛉, ψ, λ₋, Δ),
                   ∂V₋_c(𝛉, ψ, λ₋, Δ),
                   ∂V₋_d(𝛉, ψ, λ₋, Δ),
                   ∂V₋_e(𝛉, ψ, λ₋, Δ),
                   ∂V₋_f(𝛉, ψ, λ₋, Δ))
end

function ∂V₋_a(𝛉::AbstractVector, ψ::Number, λ₋::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    part1 = 1/(2*λ₋*Δ)
    part2 = (ψ/(λ₋*Δ))^(-0.5)
    part3 = 4* c*f - e^2 + 4*Δ^(-1)*c*ψ
    part4 = (ψ/(2*λ₋))*(1 - ((c - a)/((a - c)^2 + b^2)^0.5))
    return part1*part2*(part3 - part4)
end

function ∂V₋_b(𝛉::AbstractVector, ψ::Number, λ₋::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    part1 = 1/(2*λ₋*Δ)
    part2 = (ψ/(λ₋*Δ))^(-0.5)
    part3 = d*e - 2*b*f - 2*Δ^(-1)*b*ψ
    part4 = (b*ψ)/((2*λ₋)*((a - c)^2 + b^2)^0.5)
    return part1*part2*(part3 - part4)
end

function ∂V₋_c(𝛉::AbstractVector, ψ::Number, λ₋::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    part1 = 1/(2*λ₋*Δ)
    part2 = (ψ/(λ₋*Δ))^(-0.5)
    part3 = 4*a*f - d^2 + 4*Δ^(-1)*a*ψ
    part4 = (ψ/(2*λ₋))*(1 - ((a - c)/((a - c)^2 + b^2)^0.5))
    return  part1*part2*(part3 - part4)
end

function ∂V₋_d(𝛉::AbstractVector, ψ::Number, λ₋::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    part1 = (ψ/(λ₋*Δ))^(0.5)
    part2 = (b*e - 2*c*d)/(2*ψ)
    return part1*part2
end

function ∂V₋_e(𝛉::AbstractVector, ψ::Number, λ₋::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    part1 = (ψ/(λ₋*Δ))^(0.5)
    part2 = (b*d - 2*a*e)/(2*ψ)
    return part1*part2
end

function ∂V₋_f(𝛉::AbstractVector, ψ::Number, λ₋::Number, Δ::Number)
    a, b, c, d, e, f = 𝛉
    part1 = -(2*λ₋)^-1
    part2 = (ψ/(λ₋*Δ))^(-0.5)
    return part1*part2
end
