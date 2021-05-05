"""
```
    Observations
```
    A one-dimensional array of 2D points wrapped in a tuple called `data`.

    The reason I wrap the points in a tuple is that I intend to refactor `Observations` into a separate package
    in the future, and different elements of the tuple may represent corresponding points between images, for example. 
"""
struct Observations{N,T} <: AbstractObservations
    data::NTuple{N,T}
end

# Convert an N x 2 matrix to a one-dimensional array of 2D points (i.e. each point is a row of the N x 2 matrix) and wrap in a tuple.
Observations(data::AbstractMatrix; colwise = true) = Observations(tuple(svectors(transpose(data), Val{2}())))

"""
```
    UncertainObservations
```
    A one-dimensional array of 2D points wrapped in a tuple and a one-dimensional array of corresponding 2D covariance matrices (also wrapped in a tuple). 

    The reason I wrap the points in a tuple is that I intend to refactor `UncertainObservations` into a separate package
    in the future, and different elements of the tuple may represent corresponding points between images, for example. 
"""
struct UncertainObservations{N,T₁, T₂} <: AbstractUncertainObservations
    data::NTuple{N,T₁}
    covariance_matrices::NTuple{N,T₂}
end

