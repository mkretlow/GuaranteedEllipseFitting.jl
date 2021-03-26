function (::DirectEllipseFit)(observations::AbstractObservations)
    @unpack data = observations
    # # Convert the N x 2 matrix to a one-dimensional array of 2D points. 
    # ℳ = svectors(transpose(data[1]), Val{2}())
    ℳ = data[1]
    N = length(ℳ)
    𝐃₁ = zeros(N, 3)
    𝐃₂ = zeros(N, 3)
    for n = 1:N
        𝐦 = ℳ[n]
        # Quadratic part of the design matrix.
        𝐃₁[n, 1] = 𝐦[1]^2
        𝐃₁[n, 2] = 𝐦[1]*𝐦[2]
        𝐃₁[n, 3] = 𝐦[2]^2
        # Linear part of the design matrix.
        𝐃₂[n, 1] = 𝐦[1]
        𝐃₂[n, 2] = 𝐦[2]
        𝐃₂[n, 3] = 1.0
    end
    # Quadratic part of the scatter matrix.
    𝐒₁ = 𝐃₁' * 𝐃₁ 
    # Combined part of the scatter matrix.
    𝐒₂ = 𝐃₁' * 𝐃₂ 
    # Linear part of the scatter matrix.
    𝐒₃ = 𝐃₂' * 𝐃₂ 
    # For getting a₂ from a₁.
    𝐓 = -inv(𝐒₃) * 𝐒₂' 
    # Reduce scatter matrix.
    𝐌 = 𝐒₁ + 𝐒₂ * 𝐓
    # Premultiply by inv(C₁).
    𝐌₀ = vcat(𝐌[3,:]' / 2, -𝐌[2,:]', 𝐌[1,:]' / 2) 
    # Solve eigensystem.
    F = eigen(𝐌₀) 
    evec = F.vectors
    evalue = F.values
    # Evaluate a'Ca.
    cond = 4 * evec[1,:] .* evec[3,:] - evec[2,:].^2
    index = findfirst(cond .> 0)
    # Ellipse coefficients.
    𝛉₁ = SVector(evec[:, index]...)
    𝛉 = vcat(𝛉₁, 𝐓 * 𝛉₁)
    return 𝛉 / norm(𝛉)
end