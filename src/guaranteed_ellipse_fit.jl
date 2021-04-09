function (method::GuaranteedEllipseFit)(observations::AbstractObservations)
    @unpack parametrisation = method  
    @unpack optimisation_scheme = method
    @unpack seed = optimisation_scheme

    𝛈₀ = to_latent_parameters(parametrisation, seed.𝛉)
    optimisation_scheme = @set optimisation_scheme.seed = ManualEstimation(𝛈₀)
    N = length(observations.data[1])
    jacobian_matrix = JacobianMatrix(observations, zeros(N,5)) 
    𝛈₁, cost = optimisation_scheme(observations, jacobian_matrix)    
    # TODO Remove parametrisation option and just always use "SecondLatentEllipseParametrisation" 
    # The lm method below assumes "SecondLatentEllipseParametrisation"
    𝛉 = from_latent_parameters(parametrisation, 𝛈₁)   
    return 𝛉 / norm(𝛉)
end

function (lm::LevenbergMarquardt)(observations::AbstractObservations, jacobian_matrix::JacobianMatrix)
    @unpack data = observations
    @unpack γ, γ_inc, γ_dec = lm
    @unpack λ, λ_lower_bound = lm
    @unpack tol_𝛉, tol_cost, tol_Δ, tol_∇ = lm
    @unpack max_iter = lm
    @unpack seed = lm

    # minimum allowable magnitude of conic determinant to prevent ellipse from
    # convering on a degenarate parabola (e.g. two parallel lines).
    tol_detD = 1e-5
    # barrier tolerance tp prevent ellipse from convering on parabola.
    tol_barrier = 15.5

    𝛈ₖ = seed.𝛉

    keep_going = true
    theta_updated = false
    max_func_eval = 100 * length(𝛈ₖ)
    func_eval = 0
    k = 1

    ∂𝐫 = jacobian_matrix
    𝛉ₖ = from_latent_parameters(SecondLatentEllipseParametrisation(), 𝛈ₖ)
    Δₖ = ones(length(𝛉ₖ))

    N = length(first(data))
    𝐫 = vector_valued_objective(observations, 𝛉ₖ)
    # Sum of squared residuals: 𝐫'*𝐫
    costₖ = sum(abs2, 𝐫)

    𝐈 = UniformScaling(1)
    𝐅₀ = SA_F64[0 0 2; 0 -1 0; 2 0 0]
    𝐅 = SA_F64[1 0; 0 0] ⊗ 𝐅₀

    was_updated = false
    while keep_going && k < max_iter

        𝛉ₖ = from_latent_parameters(SecondLatentEllipseParametrisation(), 𝛈ₖ)
        𝐫 = vector_valued_objective(observations, 𝛉ₖ)
        𝐉 = ∂𝐫(𝛈ₖ)
        𝐇 = 𝐉'*𝐉       

        𝐏ₜ = UniformScaling(1) - (𝛉ₖ*𝛉ₖ')/norm(𝛉ₖ)^2
        ∂π = norm(𝛉ₖ)^-1 * 𝐏ₜ *  ∂𝛋(𝛈ₖ)
        𝐖 = (∂π'*∂π) 

        # Compute two potential updates based on different weightings of 𝐖.
        λ₁ = max(λ,  λ_lower_bound)
        Δ₁ = -(𝐇 + λ₁*𝐖) \ (𝐉'*𝐫)
        𝛈₁ = 𝛈ₖ + Δ₁        

        λ₂ = max(λ / γ,  λ_lower_bound)
        Δ₂ = -(𝐇 + λ₂*𝐖) \ (𝐉'*𝐫)
        𝛈₂ = 𝛈ₖ + Δ₂

        𝛉₁ = from_latent_parameters(SecondLatentEllipseParametrisation(), 𝛈₁)
        𝛉₂ = from_latent_parameters(SecondLatentEllipseParametrisation(), 𝛈₂)
        

        # Compute new residuals and costs based on these updates.
        𝐫₁ = vector_valued_objective(observations, 𝛉₁)
        cost₁ = sum(abs2, 𝐫₁)
        𝐫₂ = vector_valued_objective(observations, 𝛉₂)
        cost₂ = sum(abs2, 𝐫₂)

        if cost₁ >= costₖ && cost₂ >= costₖ
            # Neither potential update reduced the cost.
            was_updated = false
            𝛉ₖ₊₁ = 𝛉ₖ
            𝛈ₖ₊₁ = 𝛈ₖ
            Δₖ₊₁ = Δₖ
            costₖ₊₁ = costₖ
            λ = λ * γ # In the next iteration add more of the identity matrix.
            func_eval = func_eval + 1
        elseif cost₂ < costₖ
            # Update (2) reduced the cost function.
            was_updated = true
            𝛈ₖ₊₁ = 𝛈₂
            𝛉ₖ₊₁ = 𝛉₂
            Δₖ₊₁ = Δ₂
            costₖ₊₁ = cost₂
            λ = λ / γ  # In the next iteration add less of the identity matrix.
        else
            # Update (1) reduced the cost function.
            was_updated = true
            𝛈ₖ₊₁ = 𝛈₁
            𝛉ₖ₊₁ = 𝛉₁
            Δₖ₊₁ = Δ₁
            costₖ₊₁ = cost₁
            λ = λ  # Keep the same damping for the next iteration.
        end

        barrier = (𝛉ₖ₊₁'*𝐈*𝛉ₖ₊₁)/(𝛉ₖ₊₁'*𝐅*𝛉ₖ₊₁)
        D = SA_F64[𝛉ₖ₊₁[1] 𝛉ₖ₊₁[2]/2 𝛉ₖ₊₁[4]/2 ;
                   𝛉ₖ₊₁[2]/2 𝛉ₖ₊₁[3] 𝛉ₖ₊₁[5]/2 ;
                   𝛉ₖ₊₁[4]/2 𝛉ₖ₊₁[5]/2 𝛉ₖ₊₁[6]]
    
        detD = det(D)

        # Since 𝛈 is a projective entity this converge criterion will have
        # to change to take into account the scale/sign ambiguity.
        if min(norm(𝛈ₖ₊₁ - 𝛈ₖ), norm(𝛈ₖ₊₁ + 𝛈ₖ)) < tol_𝛉 && was_updated
            @info "Breaking because of tolerance."
            keep_going = false
        elseif abs(costₖ₊₁ - costₖ) < tol_cost && was_updated
            @info "Breaking because of cost."
            keep_going = false
        elseif norm(Δₖ₊₁) < tol_Δ
            @info "Breaking because of update norm."
            keep_going = false
        elseif  norm(𝐉'*𝐫, Inf) < tol_∇
            @info "Breaking because of gradient norm."
            keep_going = false
        elseif func_eval > max_func_eval
            @info "Breaking because maximum func evaluations reached."
            keep_going = false
        elseif log(barrier) > tol_barrier || abs(detD) < tol_detD
            @info "Breaking because approaching degenerate ellipse."
            keep_going = false
        end

        𝛈ₖ = 𝛈ₖ₊₁
        Δₖ = Δₖ₊₁
        costₖ = costₖ₊₁
        k = was_updated ? k + 1 : k
    end
    # TODO create optimisation output struct and return that instead
    return 𝛈ₖ, costₖ #, Δₖ, k, func_eval
end


function vector_valued_objective(observations::AbstractObservations, 𝛉::AbstractVector)
    @unpack data = observations
    ℳ = data[1]
    N = length(ℳ)
    𝐫 = zeros(N)
    for n = 1:N
        𝐦 = ℳ[n]
        𝐮ₙ = SA_F64[𝐦[1]^2, 𝐦[1]*𝐦[2], 𝐦[2]^2, 𝐦[1], 𝐦[2], 1]
        ∂𝐮ₙ = SA_F64[2*𝐦[1]  𝐦[2]  0  1  0  0; 0 𝐦[1] 2*𝐦[2] 0 1 0]'
        𝐀ₙ = 𝐮ₙ * 𝐮ₙ'
        𝚲ₙ =  SA_F64[1 0 ; 0 1]
        𝐁ₙ = ∂𝐮ₙ * 𝚲ₙ * ∂𝐮ₙ'
        𝐫[n] = sqrt(abs((𝛉' * 𝐀ₙ * 𝛉)/(𝛉' * 𝐁ₙ * 𝛉)))
    end
    return 𝐫
end


function to_latent_parameters(::FirstLatentEllipseParametrisation, 𝛉::AbstractVector)
    p = 𝛉[2] / (2*𝛉[1])
    q = ((𝛉[3] / 𝛉[1]) - (𝛉[2]/(2*𝛉[1]))^2)^(-0.5)
    r = (𝛉[4] / 𝛉[1])
    s = (𝛉[5] / 𝛉[1])
    t = (𝛉[6] / 𝛉[1])
    return SVector(p, q, r, s, t)
end

function from_latent_parameters(::FirstLatentEllipseParametrisation, 𝛋::AbstractVector)
    p, q, r, s, t = 𝛋
    a = 1
    b = 2*p
    c = p^2 + q^(-2)
    d = r
    e = s
    f = t
    𝛉 = SVector(a, b, c, d, e ,f)
    return 𝛉 / norm(𝛉)
end

function to_latent_parameters(::SecondLatentEllipseParametrisation, 𝛉::AbstractVector)
    p = 𝛉[2] / (2*𝛉[1])
    q = ((𝛉[3] / 𝛉[1]) - (𝛉[2]/(2*𝛉[1]))^2)^(0.5)
    r = (𝛉[4] / 𝛉[1])
    s = (𝛉[5] / 𝛉[1])
    t = (𝛉[6] / 𝛉[1])
    return SVector(p, q, r, s, t)
end

function from_latent_parameters(::SecondLatentEllipseParametrisation, 𝛋::AbstractVector)
    p, q, r, s, t = 𝛋
    a = 1
    b = 2*p
    c = p^2 + q^(2)
    d = r
    e = s
    f = t
    𝛉 = SVector(a, b, c, d, e ,f)
    return 𝛉 / norm(𝛉)
end

function ∂𝛋(𝛈)
    return SA_F64[0      0       0 0 0 ;
                  2      0       0 0 0 ;
                  2*𝛈[1] 2*𝛈[2]  0 0 0 ;
                  0      0       1 0 0 ;
                  0      0       0 1 0 ;
                  0      0       0 0 1]
end

function (jacobian_matrix::JacobianMatrix)(𝛈::AbstractVector)
    @unpack observations = jacobian_matrix
    @unpack data = observations
    ℳ = data[1]
    N = length(ℳ)
   
    𝛉 = from_latent_parameters(SecondLatentEllipseParametrisation(), 𝛈)
   
    𝐏ₜ = UniformScaling(1) - (𝛉*𝛉')/norm(𝛉)^2
    ∂π = norm(𝛉)^-1 * 𝐏ₜ *  ∂𝛋(𝛈)
    # TODO overwrite the pre-allocated array instead
    ∂𝐫′ = zeros(N, 5)
    for n = 1:N
        𝐦 = ℳ[n]
        𝐮ₙ = SA_F64[𝐦[1]^2, 𝐦[1]*𝐦[2], 𝐦[2]^2, 𝐦[1], 𝐦[2], 1]
        ∂𝐮ₙ = SA_F64[2*𝐦[1]  𝐦[2]  0  1  0  0; 0 𝐦[1] 2*𝐦[2] 0 1 0]'
        𝐀ₙ = 𝐮ₙ * 𝐮ₙ'
        𝚲ₙ =  SA_F64[1 0 ; 0 1]
        𝐁ₙ = ∂𝐮ₙ * 𝚲ₙ * ∂𝐮ₙ'
        𝐌ₙ = 𝐀ₙ / (𝛉' * 𝐁ₙ * 𝛉)
        𝐗ₙ =  𝐌ₙ - 𝐁ₙ * ((𝛉' * 𝐀ₙ * 𝛉)/ (𝛉' * 𝐁ₙ * 𝛉)^2)
        ∂𝐫 = (𝐗ₙ*𝛉 / sqrt(abs((𝛉' * 𝐀ₙ * 𝛉) / (𝛉' * 𝐁ₙ * 𝛉)) + eps()))'
        ∂𝐫′[n,:] = ∂𝐫 * ∂π 
     end
    return ∂𝐫′
end