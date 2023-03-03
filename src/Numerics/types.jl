# Variable dimensions
abstract type VarDim{S} end

abstract type Var{name,S<:VarDim,T,units,domain} end

abstract type AbstractDiscretization{Q,N} <: DenseArray{Q,N} end

abstract type GridSpec end
struct Edges <: GridSpec end
struct Cells <: GridSpec end

abstract type Geometry end
struct UnitVolume <: Geometry end
