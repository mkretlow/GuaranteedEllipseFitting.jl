function fit_ellipse(observations::AbstractObservations, method::AbstractFittingMethod)
    normalization_method = IsotropicScalingTranslation()
    normalize = NormalizeDataContext(observations, normalization_method)
    # Estimate ellipse in a data-driven normalised coordinate system.
    𝛉′ =  method(normalize(observations)) 
    # Transform parameters back to the original coordinate system.   
    𝛉 = normalize(FromNormalizedSpace(), 𝛉′)
    return 𝛉 / norm(𝛉)
end

function fit_ellipse(data::AbstractArray, method::AbstractFittingMethod)
    # Convert the N x 2 matrix to a one-dimensional array of 2D points. 
    ℳ = svectors(transpose(data), Val{2}())
    fit_ellipse(Observations(tuple(ℳ)), method)
end