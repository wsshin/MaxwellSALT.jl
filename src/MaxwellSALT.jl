module MaxwellSALT

using Reexport
@reexport using SALTBase, MaxwellFDM
export DirectMaxwellData

const Float = typeof(0.0)  # use Float = Float128 for quadruple precision in the future
const CFloat = Complex{Float}

const AbsVec = AbstractVector
const AbsMat = AbstractMatrix

const AbsVecComplex = AbsVec{CFloat}
const AbsMatComplex = AbsMat{CFloat}

const AbsVecNumber = AbsVec{<:Number}
const AbsMatNumber = AbsMat{<:Number}

mutable struct DirectMaxwellData{MC<:AbsMatComplex} <: LinearSolverData
    A::MC  # Maxwell operator
    CC::MC  # double curl operator ∇ × μ⁻ ¹∇ ×
    function DirectMaxwellData{MC}(A::AbsMatNumber, CC::AbsMatNumber) where {MC<:AbsMatComplex}
        mA, nA = size(A)
        mA==nA || throw(ArgumentError("A must be square."))

        mCC, nCC = size(CC)
        mCC==nCC || throw(ArgumentError("CC must be square."))

        mA==mCC || throw(ArgumentError("A and CC must have the same size."))

        new(A, CC)
    end
end
DirectMaxwellData(CC::MC) where {MC<:AbsMatComplex} = DirectMaxwellData{MC}(similar(CC), CC)  # CC is shared between modes, but A is not

Base.similar(lsd::DirectMaxwellData) = DirectMaxwellData(lsd.CC)
Base.size(lsd::DirectMaxwellData) = size(lsd.CC)

function SALTBase.init_lsd!(lsd::DirectMaxwellData, ω::Number, ε::AbsVecNumber)
    mA = size(lsd.A, 1)
    nε = length(ε)
    mA==nε || throw(ArgumentError("size(lsd.A,1) = $mA and length(ε) = $nε must be the same."))

    lsd.A .= lsd.CC  # initialize; works for sparse matrices with same nonzero entry pattern
    # info("‖CC‖₁ = $(norm(CC,1)), ω = $ω, ‖ε‖ = $(norm(ε))")
    for i = 1:nε
        lsd.A[i,i] -= ω^2 * ε[i]
    end

    return nothing
end

# Later, we can store the factorization of A in DirectMaxwellData and reuse it.
function SALTBase.linsolve!(x::AbsVecComplex, lsd::DirectMaxwellData, b::AbsVecComplex)
    x .= lsd.A \ b
end

function SALTBase.linsolve_transpose!(x::AbsVecComplex, lsd::DirectMaxwellData, b::AbsVecComplex)
    x .= lsd.A.' \ b
end

function SALTBase.linapply!(b::AbsVecComplex, lsd::DirectMaxwellData, x::AbsVecComplex)
    b .= lsd.A * x
end

end # module
