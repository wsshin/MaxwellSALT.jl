export DirectMaxwellData

mutable struct DirectMaxwellData{MC<:AbsMatComplex,F<:FactComplex} <: LinearSolverData
    CC::MC  # double curl operator ∇ × μ⁻ ¹∇ ×
    A::MC  # Maxwell operator
    fact::F
    function DirectMaxwellData{MC,F}(CC::AbsMatNumber, A::AbsMatNumber) where {MC<:AbsMatComplex,F<:FactComplex}
        mCC, nCC = size(CC)
        mCC==nCC || throw(ArgumentError("CC must be square."))

        mA, nA = size(A)
        mA==nA || throw(ArgumentError("A must be square."))

        mA==mCC || throw(ArgumentError("A and CC must have the same size."))

        new(CC, A)
    end
end

# Below, CC is shared between modes, but A is not
DirectMaxwellData(CC::SparseMatComplex) = DirectMaxwellData{SparseMatComplex,SparseLUComplex}(CC, similar(CC))
DirectMaxwellData(CC::MatComplex) = DirectMaxwellData{MatComplex,DenseLUComplex}(CC, similar(CC))
DirectMaxwellData(CC::MC) where {MC<:AbsMatComplex} = DirectMaxwellData{MC,FactComplex}(CC, similar(CC))

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

    lsd.fact = lufact(lsd.A)

    return nothing
end

# Later, we can store the factorization of A in DirectMaxwellData and reuse it.
function SALTBase.linsolve!(x::AbsVecComplex, lsd::DirectMaxwellData, b::AbsVecComplex)
    A_ldiv_B!(x, lsd.fact, b)
end

function SALTBase.linsolve_transpose!(x::AbsVecComplex, lsd::DirectMaxwellData, b::AbsVecComplex)
    At_ldiv_B!(x, lsd.fact, b)
end

function SALTBase.linapply!(b::AbsVecComplex, lsd::DirectMaxwellData, x::AbsVecComplex)
    b .= lsd.A * x
end
