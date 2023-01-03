using CryoGrid
using Dates
using DimensionalData
using Test

@testset "Diagnostics" begin
    gridvals = (0.0:0.01:10.0)u"m"
    timestamps = DateTime(2000):Day(1):DateTime(2001,12,31)
    grid = Grid(gridvals)
    @testset "Active layer thickness" begin
        x = cells(grid)
        t = convert_t.(timestamps)u"s"
        T₀ = -1.0u"°C"
        A = 2.0u"°C"
        period = uconvert(u"s", 1.0u"yr")
        α = 1e-5u"m^2/s" # diffusivity
        # damped temperature signal, neglecting depth dependent phase shift
        T(z,t) = ustrip(T₀) + ustrip(A)*exp(-z*sqrt(π/(α*period)))*cos(2π/period*t)
        T_obs = DimArray(reduce(hcat, map(xᵢ -> T.(xᵢ, t), x))u"°C", (Ti(timestamps),Z(collect(x))))
        alt = Diagnostics.active_layer_thickness(T_obs)
        true_alt = (-log(1/2)/ustrip(sqrt(π/(α*period))))u"m"
        # check that difference is less than 1 cm
        @test all(abs.(alt .- true_alt) .< 1u"cm")
    end
    @testset "Zero annual amplitude" begin
        x = cells(grid)
        t = convert_t.(timestamps)u"s"
        T₀ = -1.0u"°C"
        A = 2.0u"°C"
        period = uconvert(u"s", 1.0u"yr")
        α = 1e-7u"m^2/s" # diffusivity
        threshold = 0.1u"K"
        # damped temperature signal, neglecting depth dependent phase shift
        T(z,t) = ustrip(T₀) + ustrip(A)*exp(-z*sqrt(π/(α*period)))*cos(2π/period*t)
        T_obs = DimArray(reduce(hcat, map(xᵢ -> T.(xᵢ, t), x))u"°C", (Ti(timestamps),Z(collect(x))))
        zaa = Diagnostics.zero_annual_amplitude(T_obs; threshold)
        # divide by 4 since peak-to-peak amplitude is 2*A
        true_zaa = (-log(ustrip(threshold)/4)/ustrip(sqrt(π/(α*period))))u"m"
        # check that difference is less than 1 cm
        @test all(abs.(zaa .- true_zaa) .< 1u"cm")
    end
end
