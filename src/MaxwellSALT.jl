module MaxwellSALT

using Reexport
@reexport using SALTBase, MaxwellFDM

const Float = typeof(0.0)  # use Float = Float128 for quadruple precision in the future
const CFloat = Complex{Float}

const AbsVec = AbstractVector
const AbsMat = AbstractMatrix

const AbsVecComplex = AbsVec{CFloat}
const AbsMatComplex = AbsMat{CFloat}

const AbsVecNumber = AbsVec{<:Number}
const AbsMatNumber = AbsMat{<:Number}

const MatComplex = Matrix{CFloat}
const SparseMatComplex = SparseMatrixCSC{CFloat,Int64}

const FactComplex = Factorization{CFloat}

const SparseLUComplex = Base.SparseArrays.UMFPACK.UmfpackLU{CFloat,Int64}
const DenseLUComplex = Base.LinAlg.LU{CFloat,MatComplex}

include("direct.jl")

end # module
