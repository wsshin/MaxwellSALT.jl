module MaxwellSALT

using Reexport
@reexport using SALTBase, MaxwellFDM
using LinearAlgebra
using SparseArrays: SparseMatrixCSC
using SuiteSparse.UMFPACK: UmfpackLU

const Float = typeof(0.0)  # use Float = Float128 for quadruple precision in the future
const CFloat = Complex{Float}

const AbsVec = AbstractVector
const AbsMat = AbstractMatrix

const AbsVecComplex = AbsVec{CFloat}
const AbsMatComplex = AbsMat{CFloat}

const AbsVecReal = AbsVec{<:Real}
const AbsVecNumber = AbsVec{<:Number}
const AbsMatNumber = AbsMat{<:Number}

const MatComplex = Matrix{CFloat}
const SparseMatComplex = SparseMatrixCSC{CFloat,Int64}

const FactComplex = Factorization{CFloat}

const SparseLUComplex = UmfpackLU{CFloat,Int}
const DenseLUComplex = LU{CFloat,MatComplex}

include("direct.jl")
include("gainobj.jl")

end # module
