using CryoGrid.Numerics
using CryoGrid.Numerics: flux!, divergence!, nonlineardiffusion!
using LinearAlgebra
using Test

function test_flux!_and_divergence!()
	x = 0.0:0.25:1.0
	xc = (x[1:end-1] .+ x[2:end])./2
	Δx = x[2:end] .- x[1:end-1]
	Δxc = xc[2:end] .- xc[1:end-1]
	y = [0.0,1.0,0.0]
	k = ones(length(x))
	jy = zeros(length(x))
	flux!(jy, y, Δxc, k)
	@test jy[1] == jy[end] == zero(eltype(jy))
	@test jy[2] ≈ -1/0.25
	@test jy[3] ≈ 1/0.25
	div = zeros(length(y))
	divergence!(div, jy, Δx)
	@test allequal(div, [1/0.25^2,-2/0.25^2,1/0.25^2])
end

function test_nonlineardiffusion!()
	f(x) = (1/6)x^3 # function to differntiate
	df(x) = (1/2)x^2 # analytical 1st derivative
	d2f(x) = x # analytical 2nd derivative
	k(x) = sin(x) # non-linear diffusion function
	dk(x) = cos(x) # analytical 1st derivative of k
	d2kf(x) = k(x)d2f(x) + dk(x)df(x) # by product rule
	x = exp.(0.0:0.001:0.05)
	xc = (x[1:end-1].+x[2:end])/2
	y = f.(xc)
	kₓ = k.(x)
	d2y = d2kf.(xc) # analytical solution
	Δx = x[2:end] .- x[1:end-1]
	jy = zeros(length(x))
	# manually set boundary fluxes
	jy[1] = -kₓ[1]*df(x[1])
	jy[end] = -kₓ[end]*df(x[end])
	Δxc = xc[2:end] .- xc[1:end-1]
	div = zeros(length(y))
	nonlineardiffusion!(div, jy, y, Δxc, kₓ, Δx)
	@test allfinite(div)
	@test allequal(div, d2y, atol=0.01)
end

@testset "Math" begin
	@testset "flux! and divergence!" begin
		test_flux!_and_divergence!()
	end
	@testset "Non-linear diffusion" begin
		test_nonlineardiffusion!()
	end
	@testset "TDMA solve" begin
		A = Tridiagonal(ones(4), 2*ones(5), ones(4))
		b = ones(5)
		x_true = A \ b # solve using built-in ldiv
		x_test = similar(x_true)
		Numerics.tdma_solve!(x_test, diag(A,-1), diag(A), diag(A,1), b)
		@test all(x_true .≈ x_test)
	end
end
