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

mutable struct DirectMaxwellData{MC<:AbsMatComplex,F<:FactComplex,VC<:AbsVecComplex} <: LinearSolverData
    CC::MC  # double curl operator ∇ × μ⁻ ¹∇ ×
    A::MC  # Maxwell operator
    Mc::MC  # operator interpolating Ex, Ey, Ez at grid cell corners
    Ml::MC  # operator interpolating voxel-corner Ex, Ey, Ez at Ez, Ex, Ey (i.e., left component) locations
    Mr::MC  # operator interpolating voxel-corner Ex, Ey, Ez at Ey, Ez, Ex (i.e., right component) locations
    vtemp::VC  # temporary storage for complex vector
    wtemp::VC  # temporary storage for complex vector
    fact::F  # matrix factors of A
    function DirectMaxwellData{MC,F,VC}(CC::AbsMatNumber, A::AbsMatNumber,
                                        Mc::AbsMatNumber, Ml::AbsMatNumber, Mr::AbsMatNumber,
                                        vtemp::AbsVecNumber, wtemp::AbsVecNumber) where {MC<:AbsMatComplex,F<:FactComplex,VC<:AbsVecComplex}
        mCC, nCC = size(CC)
        mCC==nCC || throw(ArgumentError("CC must be square: size(CC) = $((mCC,nCC))."))

        mA, nA = size(A)
        (mA==mCC && nA==nCC) || throw(ArgumentError("A and CC must be the same size: size(A) = $((mA,nA)), size(CC) = $((mCC,nCC))."))

        mMc, nMc = size(Mc)
        (mMc==mCC && nMc==nCC) || throw(ArgumentError("Mc and CC must be the same size: size(Mc) = $((mMc,nMc)), size(CC) = $((mCC,nCC))."))

        mMl, nMl = size(Ml)
        (mMl==mCC && nMl==nCC) || throw(ArgumentError("Ml and CC must be the same size: size(Ml) = $((mMl,nMl)), size(CC) = $((mCC,nCC))."))

        mMr, nMr = size(Mr)
        (mMr==mCC && nMr==nCC) || throw(ArgumentError("Mr and CC must be the same size: size(Mr) = $((mMr,nMr)), size(CC) = $((mCC,nCC))."))

        new(CC, A, Mc, Ml, Mr, vtemp, wtemp)
    end
end

# Below, CC is shared between modes, but A is not (it is overwritten for new modes).
# Therefore we need to keep a separate copy of A.
DirectMaxwellData(CC::SparseMatComplex, Mc::SparseMatComplex, Ml::SparseMatComplex, Mr::SparseMatComplex) =
    (n = size(CC,2); DirectMaxwellData{SparseMatComplex,SparseLUComplex,VecComplex}(CC, similar(CC), Mc, Ml, Mr, VecComplex(undef,n), VecComplex(undef,n)))
DirectMaxwellData(CC::MatComplex, Mc::MatComplex, Ml::MatComplex, Mr::MatComplex) =
    (n = size(CC,2); DirectMaxwellData{MatComplex,DenseLUComplex,VecComplex}(CC, similar(CC), Mc, Ml, Mr, VecComplex(undef,n), VecComplex(undef,n)))
# DirectMaxwellData(CC::MC, Mc::MC, Ml::MC, Mr::MC) where {MC<:AbsMatComplex} =
#     DirectMaxwellData{MC,FactComplex}(CC, similar(CC), Mc, Ml, Mr)

Base.similar(lsd::DirectMaxwellData) = DirectMaxwellData(lsd.CC, lsd.Mc, lsd.Ml, lsd.Mr)
Base.size(lsd::DirectMaxwellData) = size(lsd.CC)

# Factorize the matrix A and cache it.
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

function SALTBase.linsolve!(x::AbsVecComplex, lsd::DirectMaxwellData, b::AbsVecNumber)
    t = @elapsed ldiv!(x, lsd.fact, b)
    # @info "time linsolve!: $t"

    return nothing
end

function SALTBase.linsolve_transpose!(x::AbsVecComplex, lsd::DirectMaxwellData, b::AbsVecNumber)
    t = @elapsed ldiv!(x, transpose(lsd.fact), b)
    # @info "time linsolve_transpose!: $t"

    return nothing
end

function SALTBase.linapply!(b::AbsVecComplex, lsd::DirectMaxwellData, x::AbsVecNumber)
    mul!(b, lsd.A, x)

    return nothing
end

function SALTBase.abs2_interp!(abs2ψ::AbsVecFloat, lsd::DirectMaxwellData, ψ::AbsVecNumber)
    # Calculate the amplitude squares of Ex, Ey, Ez at their Yee grid locations.
    abs2ψ .= abs2.(ψ)

    # Interpolate Ex, Ey, Ez at the Yee grid voxel corner locations.
    ψc = lsd.vtemp
    mul!(ψc, lsd.Mc, ψ)

    # Calculate the amplitude squares of Ex, Ey, Ez at the Yee grid's Ez, Ex, Ey locations.
    ψl = lsd.wtemp
    mul!(ψl, lsd.Ml, ψc)
    abs2ψ .+= abs2.(ψl)

    # Calculate the amplitude squares of Ex, Ey, Ez at the Yee grid's Ey, Ez, Ex locations.
    ψr = lsd.wtemp
    mul!(ψr, lsd.Mr, ψc)
    abs2ψ .+= abs2.(ψr)
end

# In IterativeMaxwellData, the only things I need to store would be the things that I need
# for applying the Maxwell operation or averaging operation matrix-free, such as ε, ω, and
# grid.  The functions themselves that apply the matrix-free Maxwell operation or
# matrix-free averaging operation to a column vector should be implemented in MaxwellFDM.
# linsolve!, etc that take IterativeMaxwellData simply need to call those functions with the
# arguments it stores.
