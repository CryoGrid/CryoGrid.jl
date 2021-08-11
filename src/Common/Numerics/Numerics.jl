module Numerics

import Base.==
import ExprTools

using CryoGrid.Utils

using Base: @inbounds, @propagate_inbounds
using DimensionalData: DimArray, dims
using IfElse
using Interpolations: Interpolations, Gridded, Linear, Flat, Line, interpolate, extrapolate
using IntervalSets
using LinearAlgebra
using LoopVectorization
using RuntimeGeneratedFunctions
using Unitful
using Symbolics
using SymbolicUtils

export ∇
export Var, Prognostic, Algebraic, Diagnostic, Parameter
export VarDim, OnGrid, Shape, Scalar
export varname, vartype, isprognostic, isalgebraic, isdiagnostic, isparameter, domain
export constrain, unconstrain

RuntimeGeneratedFunctions.init(@__MODULE__)

include("math.jl")
include("grid.jl")
include("variables.jl")

end
