
"""
```
    fit_ellipse(observations::AbstractObservations, f::AbstractFittingMethod)
```

Fits an ellipse to a set of scattered `observations` using a specifed ellipse fitting method `f`.

# Details
If each data point is accompanied by a 2x2 covariance matrix then `observations` is an instance of
[UncertainObservations](@ref), otherwise it is an instance of `[Observations](@ref)`

If `f` is `GuaranteedEllipseFit()` and `observations <: UncertainObservations` then 
the supplied covariance matrices are incorporated into the estimation procedure. 
If `f` is `GuaranteedEllipseFit()` and `observations <: Observations` (data point covariance matrices
are not provided), then the algorithm assumes homogeneous Gaussian noise and estimates the
level of noise from the data. 

If `f` is `DirectEllipseFit()` then data point covariance matrices are ignored even if they 
are supplied because the algroithm cannot incorporate such information. 
"""
function fit_ellipse(observations::AbstractObservations, method::AbstractFittingMethod)
    normalization_method = IsotropicScalingTranslation()
    normalize = NormalizeDataContext(observations, normalization_method)
    # Estimate ellipse in a data-driven normalised coordinate system.
    𝛉′ =  method(normalize(observations)) 
    # Transform parameters back to the original coordinate system.   
    𝛉 = normalize(FromNormalizedSpace(), 𝛉′)
    return 𝛉 / norm(𝛉)
end

function fit_ellipse(data::AbstractMatrix, method::AbstractFittingMethod)
    # Convert the N x 2 matrix to a one-dimensional array of 2D points. 
    ℳ = svectors(transpose(data), Val{2}())
    fit_ellipse(Observations(tuple(ℳ)), method)
end

function fit_ellipse(observations::AbstractObservations, method::GuaranteedEllipseFit)
    @unpack optimisation_scheme = method
    @unpack seed = optimisation_scheme
    normalization_method = IsotropicScalingTranslation()
    normalize = NormalizeDataContext(observations, normalization_method)
    # Convert observations and seed to a data-driven normalised coordinate system.
    method = @set method.optimisation_scheme.seed = ManualEstimation(normalize(ToNormalizedSpace(), seed.𝛉))
    # Apply the guaranteed ellipse fit. 
    𝛉′ =  method(normalize(observations))  
    # Transform parameters back to the original coordinate system.   
    𝛉 = normalize(FromNormalizedSpace(), 𝛉′)   
    return 𝛉 / norm(𝛉)
end