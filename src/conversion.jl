function (::GeometricToAlgebraic)(𝛏::AbstractVector)
    A, B, H, K, τ = 𝛏    

    a = cos(τ)^2/A^2 + sin(τ)^2/B^2
    b = (1/A^2 - 1/B^2)*sin(2*τ)
    c = (cos(τ)^2/B^2) + sin(τ)^2/A^2;
    d = (2*sin(τ)*(K*cos(τ)-H*sin(τ))) / B^2 - (2*cos(τ)^2*(H+K*tan(τ))) / A^2
    e = (2*cos(τ)*(H*sin(τ) - K*cos(τ))) / B^2 - (2*sin(τ)*(H*cos(τ) + K*sin(τ))) / A^2
    f = (H*cos(τ) + K*sin(τ))^2 / A^2 + (K*cos(τ) - H*sin(τ))^2 / B^2 - 1

    𝛉 = SVector(a, b, c, d, e, f)
    return 𝛉
end

function (::AlgebraicToGeometric)(𝛉::AbstractVector)
    a, b, c, d, e, f  = 𝛉   
    Δ = b^2 - 4*a*c
    λ₊ = 0.5*(a + c - (b^2 + (a - c)^2)^0.5)
    λ₋ = 0.5*(a + c + (b^2 + (a - c)^2)^0.5)

    ψ = b*d*e - a*e^2 - b^2*f + c*(4*a*f - d^2)
    V₊ = (ψ/(λ₊*Δ))^0.5
    V₋ = (ψ/(λ₋*Δ))^0.5

    # major semi-axis
    A = max(V₊, V₋)
    # minor semi-axis
    B = min(V₊, V₋)

    # determine x-coordinate of ellipse centroid
    H = (2*c*d - b*e)/(Δ)
    # determine y-coordinate of ellipse centroid
    K = (2*a*e - b*d)/(Δ)

    # angle between x-axis and major axis 
    τ = 0
    # determine tilt of ellipse in radians
    if V₊ >= V₋
        if (b == 0 && a < c)
            τ = 0
        elseif (b == 0 && a >= c)
            τ  = 0.5*π
        elseif (b < 0 && a < c)
            τ  = 0.5*acot((a - c)/b)
        elseif (b < 0 && a == c) 
            τ = π/4
        elseif (b < 0 && a > c)
            τ = 0.5*acot((a - c)/b) + π/2
        elseif (b > 0 && a < c)
            τ = 0.5*acot((a - c)/b) + π
        elseif (b > 0 && a == c)
            τ = π*(3/4)
        elseif (b > 0 && a > c)
            τ = 0.5*acot((a - c)/b) + π/2
        end
    elseif V₊ < V₋
        if (b == 0 && a < c)
            τ = π/2
        elseif (b == 0 && a >= c)
            τ = 0
        elseif (b < 0 && a < c)
            τ = 0.5*acot((a - c)/b) + π/2
        elseif (b < 0 && a == c)
            τ = π*(3/4)
        elseif (b < 0 && a > c)
            τ = 0.5*acot((a - c)/b) + π
        elseif (b > 0 && a < c)
            τ = 0.5*acot((a - c)/b) + π/2
        elseif (b > 0 && a == c)
            τ = π/4
        elseif (b > 0 && a > c)
            τ = 0.5*acot((a - c)/b)
        end
    end
    𝛏 = SVector(A, B, H, K, τ)
    return 𝛏
end