# Below, there are two ways to define the functions with the same names as those defined
# in other packages, e.g., MaxwellFDM: extending the functions defined in the other packaces,
# or calling the functions defined in the other packages by qualifying them with the package
# names.
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

# If we are extending a function that is already exported by a reexported package, I think
# the former is better, because it explicitly specifies that we are extending the
# capabilities of an already exported function.
#
# On the other hand, if we are "extending" a function that is exported by a not-reexported
# package (e.g., SALTBase), then we need to export the function anyway to make it available
# to the users.  Then, we are actually not extending the capability of some function that is
# already available to the users.  Therefore, in this case I think using the latter makes
# more sense.
#
# Currently, I am reexporting MaxwellFDM from the present package, but not SALTBase.  I
# think the users may want to solve a driven Maxwell's equations, so I reexport MaxwellFDM.
# However, the users do not need to know the internals of SALTBase, so I do not reexport
# SALTBase.  Therefore, when redefining functions already defined in MaxwellFDM and SALTBase
# below, I extend MaxwellFDM's functions (the former), but call SALTBase's functions with
# qualification (the latter).

export SALT
export get_εcold, set_gainparam!, set_initguess!, get_nonlasingsol, get_lasingsol,
       add_gainobj!, simulate!

mutable struct SALT
    m::Maxwell

    # Cold materials
    εc::AbsVecComplex

    # Gain parameters
    ω₀::Real
    γperp::Real
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
MaxwellFDM.set_Npml!(s::SALT, Npml::Tuple2{AbsVecInteger}) = set_Npml!(s.m, Npml)

MaxwellFDM.set_background!(s::SALT, matname::String, ε::MatParam) = set_background!(s.m, matname, ε)
MaxwellFDM.add_obj!(s::SALT, matname::String, ε::MatParam, shapes::Shape...) = add_obj!(s.m, matname, ε, shapes...)
MaxwellFDM.add_obj!(s::SALT, matname::String, ε::MatParam, shapes::AbsVec{<:Shape}) = add_obj!(s.m, matname, ε, shapes)

#= Getters delegating to Maxwell =#
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
set_gainparam!(s::SALT, ω₀::Real, γperp::Real) = (s.ω₀ = ω₀; s.γperp = γperp; set_freq!(s.m, -ω₀))

function get_gainprofile(s::SALT)
    if ~isdefined(s, :gp)
        g = get_grid(s.m)
        N = 3*prod(g.N)
        s.gp = GainProfile(s.ω₀, s.γperp, N)
    end

    return s.gp
end

# Later, when iterative solvers are supported, I will need to parametrize SALT as
# SALT{<:MaxwellData}, such that the users can select which class of solvers to use.
function get_maxwelldata(s::SALT)
    if ~isdefined(s, :lsd)
        CC = get_dblcurl(s)
        s.lsd = DirectMaxwellData(CC)
    end

    return s.lsd
end

set_initguess!(s::SALT, ωguess::AbsVecNumber, Ψguess::AbsMatComplex) = (s.ωguess = ωguess; s.Ψguess = Ψguess)

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

        g = get_grid(s.m)
        N = 3*prod(g.N)
        M = length(s.ωguess)

        s.nlvar = NonlasingVar(lsd, N, M)
    end

    return s.nlvar
end

function get_lasingsol(s::SALT)
    if ~isdefined(s, :lsol)
        g = get_grid(s.m)
        N = 3*prod(g.N)
        M = length(s.ωguess)

        s.lsol = LasingSol(N, M)
    end

    return s.lsol
end

function get_lasingvar(s::SALT)
    if ~isdefined(s, :lvar)
        lsd = get_maxwelldata(s)

        g = get_grid(s.m)
        N = 3*prod(g.N)
        M = length(s.ωguess)

        s.lvar = LasingVar(lsd, N, M)
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
    g = get_grid(s.m)
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
