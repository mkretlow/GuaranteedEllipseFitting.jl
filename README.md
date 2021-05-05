# GuaranteedEllipseFitting

<div class="row">
  <div class="column">
   <img src="https://github.com/zygmuntszpak/GuaranteedEllipseFitting.jl/blob/master/docs/src/images/homogeneous_noise_example.gif"/>
  </div>
</div>

A Julia implementation of the paper 

>Szpak, Z. L., Chojnacki, W., & van den Hengel, A. (2014). Guaranteed Ellipse Fitting with a Confidence Region and an Uncertainty Measure for Centre, Axes, and Orientation. Journal of Mathematical Imaging and Vision, 52(2), 173–199. doi:10.1007/s10851-014-0536-x

## Abstract
A simple and fast ellipse estimation method is presented based on optimisation of the Sampson distance serving as a measure of the quality of fit between a candidate ellipse and data points. Generation of ellipses, not just conics, as estimates is ensured through the use of a parametrisation of the set of all ellipses. Optimisation of the Sampson distance is performed with the aid of a custom variant of the Levenberg–Marquardt algorithm. The method is supplemented with a measure of uncertainty of an ellipse fit in two closely related forms. One of these concerns the uncertainty in the algebraic parameters of the fit and the other pertains to the uncertainty in the geometrically meaningful parameters of the fit such as the centre, axes, and major axis orientation. In addition, a means is provided for visualising the uncertainty of an ellipse fit in the form of planar confidence regions. For moderate noise levels, the proposed estimator produces results that are fully comparable in accuracy to those produced by the much slower maximum likelihood estimator. Due to its speed and simplicity, the method may prove useful in numerous industrial applications where a measure of reliability for geometric ellipse parameters is required.

## Usage Example

The following script reproduces a still image of the animation displayed at the start of this README.  

```julia
using CairoMakie
using Colors
using ColorSchemes
using LinearAlgebra
using GuaranteedEllipseFitting
using StaticArrays
using UnPack
using Random

# Perturb the list of 2D points ℳ according to a list of 2D covariance matrices 𝒞 
# (one covariance matrix 𝚲 for each point 𝐦 in ℳ).
function apply_noise(ℳ::AbstractVector, 𝒞::AbstractVector)
    𝒪 = similar(ℳ)
    dim = length(first(ℳ))
    for i in eachindex(ℳ)
        𝚲 = 𝒞[i]
        𝐦 = ℳ[i]
        F = cholesky(Symmetric(𝚲))
        𝐕 = F.L
        𝐫 = @SVector randn(dim)
        Δ = 𝐕*𝐫
        𝒪[i] = 𝐦 + Δ
    end
    return 𝒪
end

function homogeneous_noise_example()
    Random.seed!(5)    
    # [semi-major, semi-minor, center_x, center_y, tilt].
    𝛏 = [10, 5, 25, 25, π/4]
    # A functor to convert geometric parameters to algebraic parameters.
    # To convert algebraic to geometric, use the functor AlgebraicToGeometric(). 
    geo_to_alg = GeometricToAlgebraic()
    𝛉₀ = geo_to_alg(𝛏)
    # Sampling angles.
    α₁ = π/4
    α₂ = 1.4*π
    # Total number of points.
    N = 50
    # Refer to the docstring to learn about the `Observations` struct. 
    observations = sample_points_on_ellipse(𝛏..., N, α₁, α₂)
    # An observation has data, and the data is a one-dimensional array of 2D points that is wrapped in a tuple. 
    # We use the @unpack macro from the UnPack package to save some typing. 
    @unpack data = observations
    # pts is now a one-dimensional array of 2D points (SVector{2,Float64} from the StaticArrays package). 
    pts = data[1]    
    # Add noise to the sampled points. 
    σ² = 0.2 
    # For simplicity we shall model homogeneous and isotropic Gaussian noise. 
    # However, the guaranteed ellipse fit method is suitable for the most general case of 
    # inhomogeneous anisotropic noise. 
    𝚲 = SMatrix{2,2}(I) * σ²
    covariance_matrices = [𝚲 for n = 1:N]
    # Perturb the points so that they no longer lie on the ellipse.
    noisy_pts = apply_noise(pts, covariance_matrices)
    # UncertainObservations are characterised by the fact that each data point has an 
    # accompanying covariance matrix. 
    noisy_observations = UncertainObservations(tuple(noisy_pts), tuple(covariance_matrices))   
    
    # Direct ellipse fit estimate.
    𝛉₁ = fit_ellipse(noisy_observations, DirectEllipseFit())
    # Guaranteed ellipse fit estimate using the direct ellipse fit as a seed.
    𝛉₂ = fit_ellipse(noisy_observations, GuaranteedEllipseFit(𝛉₁))    
    
    # Evaluate algebraic ellipse equation to determine the zero level-set.
    # That is, we want to know where the equation ax^2 + bxy + cy^2 + dx + ey + f = 0 holds, 
    # where 𝛉 = [a,b,c,d,e,f]. 
    xs = LinRange(10, 40, 100)
    ys = LinRange(10, 40, 100)
    u(x,y) = SVector(x^2, x*y, y^2, x, y, 1)
    zs₀ = [dot(u(x,y), 𝛉₀) for x in xs, y in ys]
    zs₁ = [dot(u(x,y), 𝛉₁) for x in xs, y in ys]
    zs₂ = [dot(u(x,y), 𝛉₂) for x in xs, y in ys]   
    f = Figure(resolution = (800, 600))
    Axis(f[1, 1])
    
    # Plot the location of the 2D points.  
    scatter!(noisy_pts)
    # Plot the ground truth ellipse. 
    contour!(xs, ys, zs₀; levels = [0], linewidth = 2, color=:black)
    # Plot the "direct ellipse fit" result.
    contour!(xs, ys, zs₁; levels = [0], linewidth = 2, color=RGB{Float64}(230/255, 23/155, 181/255))
    # Plot the "guaranteed ellipse fit" result.
    contour!(xs, ys, zs₂; levels = [0], linewidth = 2, color=RGB{Float64}(11/255, 107/155, 105/255))  
    
    # Plot the planar confidence region.
    # We have five degrees of freedom (because the scale of 𝛉 doesn't matter), and for a p-value of
    # 0.05 the corresponding cut-off is 11.07. 
    t = 11.07
    𝚲ₐ = determine_algebraic_covariance(noisy_observations, 𝛉₂, GuaranteedEllipseFit())
    zsₐ = [(𝛉₂'*u(x,y)*u(x,y)'*𝛉₂) / (u(x,y)'*𝚲ₐ*u(x,y))  for x in xs, y in ys]
    myscheme = ColorScheme([RGBA{Float64}(0, 230/255, 28/255, 0.2)])
    # We want to fill the contour in the range [0.0,t] but the contourf
    # range requires at least 3 values, hence this "workaround". 
    contourf!(xs, ys, zsₐ; levels = [0.0, t-1, t], colormap = myscheme)
    return f
end
```

To run the above code, copy and pase it into a Julia REPL or write a script in VSCode. Note that you will first need to install all the "used" packages. 
For working with unicode, checkout https://juliamono.netlify.app/. 

