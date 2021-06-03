module MaxwellSALT

using Reexport
@reexport using MaxwellWave
using SALTBase, LinearAlgebra
using SparseArrays: SparseMatrixCSC
using SuiteSparse.UMFPACK: UmfpackLU
using AbbreviatedTypes

const SparseMatComplex = SparseMatrixCSC{CFloat,Int}
const FactComplex = Factorization{CFloat}
const SparseLUComplex = UmfpackLU{CFloat,Int}
const DenseLUComplex = LU{CFloat,MatComplex}

include("direct.jl")
include("gainobj.jl")
include("salt.jl")

end # module MaxwellSALT
