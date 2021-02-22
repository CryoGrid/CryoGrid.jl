ENV["CG_DEBUG"] = "true"
using CryoGrid
using Test

include("../testutils.jl")

@testset "∇ (fist order finite diff)" begin
	f(x) = (1/2)x^2 # function to differntiate
	df(x) = x # analytical 2nd derivative
	x = exp.(0.0:0.001:0.05)
	y = f.(x)
	dy = df.(x[1:end-1])
	δ = x[2:end] .- x[1:end-1]
	∂y = zeros(length(y)-1)
	∇(y,δ,∂y)
	@test allfinite(∂y)
	@test allequal(∂y,dy,atol=0.01)
end

@testset "∇² (second order finite diff, constant k)" begin
	f(x) = (1/6)x^3 # function to differntiate
	d2f(x) = x # analytical 2nd derivative
	x = exp.(0.0:0.001:0.05)
	y = f.(x)
	d2y = d2f.(x[2:end-1])
	k = 2.0 # constant diffusion
	δ = x[2:end] .- x[1:end-1]
	∂²y = zeros(length(y)-2)
	∇²(y,δ,k,∂²y)
	@test allfinite(∂²y)
	@test allequal(∂²y,k*d2y,atol=0.01)
end

@testset "∇² (second order finite diff, non-linear k)" begin
	f(x) = (1/6)x^3 # function to differntiate
	df(x) = (1/2)x^2 # analytical 1st derivative
	d2f(x) = x # analytical 2nd derivative
	k(x) = sin(x) # non-linear diffusion function
	dk(x) = cos(x) # analytical 1st derivative of k
	d2kf(x) = k(x)d2f(x) + dk(x)df(x) # by product rule
	x = exp.(0.0:0.001:0.05)
	y = f.(x)
	kₓ = k.((x[1:end-1] .+ x[2:end])/2)
	d2y = d2kf.(x[2:end-1]) # analytical solution
	δ = x[2:end] .- x[1:end-1]
	∂²y = zeros(length(y)-2)
	∇²(y,δ,kₓ,∂²y)
	@test allfinite(∂²y)
	@test allequal(∂²y,d2y,atol=0.01)
end
