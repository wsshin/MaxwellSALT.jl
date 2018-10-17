export DirectMaxwellData

# To do:
#
# We can create a similar type that reduces the number of factorization.  During the Newton
# iteration for solving the nonlinear equation, the matrix changes only diagonally.
# Furthermore, this diagonal change is parametrized by ω.  Therefore, we can perform
# factorization only once at a single ω, and use the idea of fast frequency sweep for other
# ω's.  Because γ(ω) has a pole, we may need to use the Padé approximation instead of the
# simple power series expansion, for which γ(ω) will generate infinitely many terms because
# γ(ω) is infinitely differentiable.
#
# A simiar idea cannot be applied to solving the lasing equation, where the diagonal change
# involves a change in the hole-burning term and is not parametrized by the scalar ω.
# However, we may construct a new parametrized matrix for each hole-burning term, by fixing
# the hole burning term and changing only ω.  The downside of this is that the coefficients
# of the series expansion of the solution needs to be calculated for every ω.  (For the
# nonlasing equation, this needed to be calculated only once for varying ω during the Newton
# iteration...  Is this right?  This is correct when both A and b are parametrized by ω, but
# in our case b is probably not, so this may not be correct.)  However, if the major
# performance improvement comes from the reduced number of factorizations rather than the
# reduced number of linear solves, this scheme should still have a huge benefit.
#
# However, all these ideas become irrelevant for 3D problems, where linear systems will be
# solved iteratively and therefore no factorization will be performed.  Hence, implementing
# the above ideas is not a priority.

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
    for i = 1:nε
        lsd.A[i,i] -= ω^2 * ε[i]
    end

    t = @elapsed lsd.fact = lu(lsd.A)
    # @info "time factorization: $t"

    return nothing
end

# Later, we can store the factorization of A in DirectMaxwellData and reuse it.
function SALTBase.linsolve!(x::AbsVecComplex, lsd::DirectMaxwellData, b::AbsVecComplex)
    t = @elapsed ldiv!(x, lsd.fact, b)
    # @info "time linsolve!: $t"

    return nothing
end

function SALTBase.linsolve_transpose!(x::AbsVecComplex, lsd::DirectMaxwellData, b::AbsVecComplex)
    t = @elapsed ldiv!(x, transpose(lsd.fact), b)
    # @info "time linsolve_transpose!: $t"

    return nothing
end

function SALTBase.linapply!(b::AbsVecComplex, lsd::DirectMaxwellData, x::AbsVecComplex)
    b .= lsd.A * x

    return nothing
end
