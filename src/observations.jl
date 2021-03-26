
struct Observations{N,T} <: AbstractObservations
    data::NTuple{N,T}
end

struct UncertainObservations{N,T₁, T₂} <: AbstractUncertainObservations
    data::NTuple{N,T₁}
    covariance_matrices::NTuple{N,T₂}
end

