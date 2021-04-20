using CryoGrid
using Test

@testset "Utils" begin
    @testset "generate_derivative" begin
        # test 1: single arg
        f(x) = x^2 + 2*x
        ∇f(x) = 2x + 2
        ∇f_gen = generate_derivative(f, :x)
        xs = -10:0.1:10
        true_vals = [∇f(x) for x in xs]
        gen_vals = [∇f_gen(x) for x in xs]
        @test all(true_vals .≈ gen_vals)
        # test 2: with extra args
        f2(x,a,b) = a*x^2 + 2*x*b
        ∇f2(x,a,b) = 2x*a + 2*b
        ∇f2_gen = generate_derivative(f2, :x)
        xs = -1:0.1:1
        as = -1:0.1:1
        bs = -1:0.1:1
        true_vals = [∇f2(x,a,b) for x in xs for a in as for b in bs]
        gen_vals = [∇f2_gen(x,a,b) for x in xs for a in as for b in bs]
        @test all(true_vals .≈ gen_vals)
        # test 3: with extra args, different order
        f3(a,x,b) = a*x^2 + 2*x*b
        ∇f3(a,x,b) = 2x*a + 2*b
        ∇f3_gen = generate_derivative(f3, :x)
        xs = -1:0.1:1
        as = -1:0.1:1
        bs = -1:0.1:1
        true_vals = [∇f3(a,x,b) for x in xs for a in as for b in bs]
        gen_vals = [∇f3_gen(a,x,b) for x in xs for a in as for b in bs]
        @test all(true_vals .≈ gen_vals)
    end
end
