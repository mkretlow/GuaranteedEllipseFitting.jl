function sample_points_on_ellipse(A::Real, B::Real, H::Real, K::Real, τ::Real, N::Integer, α₁::Real, α₂::Real)
    ℳ = [SVector(0.0,0.0) for n = 1:N]
    for (n,α) in enumerate(range(α₁, stop = α₂, length = N))
        x = H + A*cos(α)*cos(τ) - B*sin(α)*sin(τ)
        y = K + A*cos(α)*sin(τ) + B*sin(α)*cos(τ)
        𝐦 = SVector(x,y)
        ℳ[n] = 𝐦
    end
    data = tuple(ℳ)
    return Observations(data)
end