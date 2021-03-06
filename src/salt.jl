# Below, there are two ways to define the functions with the same names as those defined
# in other packages (e.g., MaxwellFDM): extending the functions defined in the other
# packaces, or calling the functions defined in the other packages by qualifying them with
# the package names.
#
# For example, for set_unitlen! below, we could define it as either
#     MaxwellFDM.set_unitlen!(s::SALT, unitlen::Real) = set_unitlen!(s.m, unitlen)
# or
#     set_unitlen!(s::SALT, unitlen::Real) = MaxwellFDM.set_unitlen!(s.m, unitlen)
#
# Both definitions works.  If we reexport MaxwellFDM from the present package, we don't even
# need to export the function, because the function name is already exported by MaxwellFDM.
# (If MaxwellFDM were not reexported, we would need to export set_unitlen! explicitly in
# this file.)
#
# Which one is the right approach?

# If we are extending a function that is already exported by the present package (e.g.,
# using reexport), I think the former is better, because it explicitly specifies that we are
# extending the capabilities of an already exported function.
#
# On the other hand, if we are "extending" a function that is exported by package (e.g.,
# SALTBase) that is not reexported by the present package, then we need to export the
# function anyway to make it available to the users.  Then, we are actually not extending
# the capability of some function that is already available to the users.  Therefore, in
# this case I think using the latter makes more sense.
#
# Currently, I am reexporting MaxwellFDM from the present package, but not SALTBase.  I
# think the users may want to solve a driven Maxwell's equations, so I reexport MaxwellFDM.
# However, the users do not need to know the internals of SALTBase, so I do not reexport
# SALTBase.  Therefore, when redefining functions already defined in MaxwellFDM and SALTBase
# below, I extend MaxwellFDM's functions (the former), but call SALTBase's functions with
# qualification (the latter).

export SALT
export set_pmlfreq!, get_εcold, set_gainparam!, set_initguess!, get_nonlasingsol,
       get_lasingsol, add_gainobj!, simulate!

mutable struct SALT
    m::Maxwell

    # Cold materials
    εc::AbsVecComplex

    # Gain parameters
    ω₀::Union{Real,AbsVecReal}
    γperp::Union{Real,AbsVecReal}
    wt::AbsVecReal
    gp::GainProfile

    # Gain objects
    gobj_vec::AbsVec{GainObject{3}}

    # Initial guess
    ωguess::AbsVecNumber
    Ψguess::AbsMatComplex

    # Solver
    lsd::DirectMaxwellData

    # Solution placeholders
    nlsol::NonlasingSol
    nlvar::NonlasingVar
    lsol::LasingSol
    lvar::LasingVar

    function SALT()
        s = new()

        s.m = Maxwell()
        s.gobj_vec = GainObject[]

        return s
    end
end


#= Setters delegating to Maxwell =#
MaxwellFDM.set_unitlen!(s::SALT, unitlen::Real) = set_unitlen!(s.m, unitlen)
MaxwellFDM.set_bounds!(s::SALT, bounds::Tuple2{AbsVecReal}) = set_bounds!(s.m, bounds)
MaxwellFDM.set_∆l!(s::SALT, ∆l::AbsVecReal) = set_∆l!(s.m, ∆l)
MaxwellFDM.set_isbloch!(s::SALT, isbloch::AbsVecBool) = set_isbloch!(s.m, isbloch)
MaxwellFDM.set_kbloch!(s::SALT, kbloch::AbsVecReal) = set_kbloch!(s.m, kbloch)
MaxwellFDM.set_Npml!(s::SALT, Npml::Tuple2{AbsVecInteger}) = set_Npml!(s.m, Npml)

MaxwellFDM.set_background!(s::SALT, matname::String, ε::MatParam) = set_background!(s.m, matname, ε)
MaxwellFDM.add_obj!(s::SALT, matname::String, ε::MatParam, shapes::Shape...) = add_obj!(s.m, matname, ε, shapes...)
MaxwellFDM.add_obj!(s::SALT, matname::String, ε::MatParam, shapes::AbsVec{<:Shape}) = add_obj!(s.m, matname, ε, shapes)

# Note that the sign of ω flips below.  This is to account for the different time-dependence
# convention between SALTBase (exp(-iωt)) and MaxwellFDM (exp(+iωt)).
#
# What is really done to obtain the physicist's Maxwell's equations from the engineer's
# Maxwell's equations is to take the complex conjugate of the entire equations.  This makes
# the time dependence from exp(+iωt) to conj{exp(+i ω_eng t)} = exp(-i conj(ω_eng) t). This
# is why we use the conjugated frequency ω_phy = conj(ω_eng) in SALT.
#
# The PML formula is embedded inside Maxwell's equations, so it is also affected by the
# conjugation.  It turns out that i appears only in the form of i ω_eng, and all the other
# PML parameters except for ω_eng are real.  Therefore, the PML formula for SALT is obtained
# by using -i instead of i and conj(ω_eng) instead of ω_eng.
#
# In SALT, we are already using conj(ω_eng) = ω_phy, so we only need to flip the sign of i
# in order to optain the PML formula for the SALT equation.  Because i appears only in the
# form of iω in the PML formula, this can be achieved by inputting -ω_phy in MaxwellFDM's
# PML formula.
set_pmlfreq!(s::SALT, ω::Number) = set_freq!(s.m, -ω)

#= Getters delegating to Maxwell =#
MaxwellFDM.get_grid(s::SALT) = get_grid(s.m)
MaxwellFDM.get_dblcurl(s::SALT) = get_dblcurl(s.m)

function get_εcold(s::SALT)
    if ~isdefined(s, :εc)
        s.εc = diag(get_εmatrix(s.m))
    end

    return s.εc
end


#= SALT's own functions =#

# Below, the minus sign is introduced because MaxwellFDM assumes the exp(+iωt) time
# dependence, whereas SALTBase assumes the exp(-iωt) time dependence.
set_gainparam!(s::SALT, ω₀::Real, γperp::Real) = (s.ω₀ = ω₀; s.γperp = γperp; return nothing)
set_gainparam!(s::SALT, ω₀::AbsVecReal, γperp::AbsVecReal, wt::AbsVecReal) = (s.ω₀ = ω₀; s.γperp = γperp; s.wt = wt; return nothing)

function get_gainprofile(s::SALT)
    if ~isdefined(s, :gp)
        g = get_grid(s)
        N = 3*prod(g.N)
        if ~isdefined(s, :wt)
            s.gp = GainProfile(s.ω₀, s.γperp, N)
        else
            s.gp = GainProfile(s.ω₀, s.γperp, N, s.wt)
        end
    end

    return s.gp
end

# Later, when iterative solvers are supported, I will need to parametrize SALT as
# SALT{<:MaxwellData}, such that the users can select which class of solvers to use.
function get_maxwelldata(s::SALT)
    if ~isdefined(s, :lsd)
        CC = get_dblcurl(s)
        Mc = get_Mc(s.m)
        Ml = get_Ml(s.m)
        Mr = get_Mr(s.m)

        s.lsd = DirectMaxwellData(CC, Mc, Ml, Mr)
    end

    return s.lsd
end

set_initguess!(s::SALT, ωguess::AbsVecNumber, Ψguess::AbsMatComplex) = (s.ωguess = ωguess; s.Ψguess = Ψguess; return nothing)

function get_nonlasingsol(s::SALT)
    if ~isdefined(s, :nlsol)
        s.nlsol = NonlasingSol(s.ωguess, s.Ψguess)
    end
    normalize!(s.nlsol)

    return s.nlsol
end

function get_nonlasingvar(s::SALT)
    if ~isdefined(s, :nlvar)
        lsd = get_maxwelldata(s)

        g = get_grid(s)
        N = 3*prod(g.N)
        M = length(s.ωguess)

        s.nlvar = NonlasingVar(lsd, N, M)
    end

    return s.nlvar
end

function get_lasingsol(s::SALT)
    if ~isdefined(s, :lsol)
        g = get_grid(s)
        N = 3*prod(g.N)
        M = length(s.ωguess)

        s.lsol = LasingSol(N, M)
    end

    return s.lsol
end

function get_lasingvar(s::SALT)
    if ~isdefined(s, :lvar)
        lsd = get_maxwelldata(s)

        g = get_grid(s)
        N = 3*prod(g.N)
        M = length(s.ωguess)
        gp = get_gainprofile(s)

        s.lvar = LasingVar(lsd, N, length(gp), M)
    end

    return s.lvar
end

#= Gain objects =#
add_gainobj!(s::SALT, shapes::Shape...) = add_gainobj!(s, d::Real->d, shapes...)
add_gainobj!(s::SALT, D₀fun::Function, shapes::Shape...) = add_gainobj!(s, D₀fun, [shapes...])
add_gainobj!(s::SALT, shapes::AbsVec{<:Shape}) = add_gainobj!(s, d::Real->d, shapes)

function add_gainobj!(s::SALT, D₀fun::Function, shapes::AbsVec{<:Shape})
    for shape = shapes  # shapes is tuple
        gobj = GainObject(shape, D₀fun)
        push!(s.gobj_vec, gobj)
    end

    return nothing
end

function simulate!(s::SALT,
                   dvec::AbsVecReal;  # trajectory of pump strength parameter to follow
                   outωaψ::NTuple{3,Bool}=(true,true,false),  # true to output ω, a, ψ
                   doutvec::AbsVecReal=dvec,  # output results when dvec[i] ∈ doutvec
                   τr_newton::Real=TR_NEWTON,  # relative tolerance for Newton method to solve nonlasing equation
                   τa_newton::Real=TA_NEWTON,  # absolute tolerance for Newton method to solve nonlasing equation
                   maxit_newton::Integer=MAXIT_NEWTON,  # maximum number of Newton iteration steps
                   m_anderson::Integer=M_ANDERSON,  # number of basis vectors to use in Anderson acceleration
                   τr_anderson::Real=TR_ANDERSON,  # relative tolerance; consider using Base.rtoldefault(Float)
                   τa_anderson::Real=TA_ANDERSON,  # absolute tolerance
                   maxit_anderson::Integer=typemax(Int),  # maximum number of Anderson iteration steps
                   verbose::Bool=true)
    g = get_grid(s)
    setD₀!(gp::GainProfile, d::Real) = assign_pumpstr!(get_gainprofile(s).D₀, s.gobj_vec, d, g.N, g.l)

    # Return, as functions of d, the followings:
    # eigenfrequencies, modal amplitudes, normalized mode profiles, numbers of AA steps, times taken for convergence of AA
    ωout, aout, ψout, nAA, tAA =
        SALTBase.simulate!(get_lasingsol(s),
                  get_lasingvar(s),
                  get_nonlasingsol(s),
                  get_nonlasingvar(s),
                  get_gainprofile(s),
                  get_εcold(s),
                  dvec,
                  setD₀!,
                  outωaψ=outωaψ,
                  doutvec=doutvec,
                  τr_newton=τr_newton,
                  τa_newton=τa_newton,
                  maxit_newton=maxit_newton,
                  m_anderson=m_anderson,
                  τr_anderson=τr_anderson,
                  τa_anderson=τa_anderson,
                  maxit_anderson=maxit_anderson,
                  verbose=verbose)

    return ωout, aout, ψout, nAA, tAA
end
