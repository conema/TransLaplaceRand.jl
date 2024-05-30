struct NewtonBracket <: Roots.AbstractBracketingMethod end
Roots.fn_argout(::NewtonBracket) = 2

struct NewtonBracketState{T,S} <: AbstractUnivariateZeroState{T,S}
    xn0::T
    xn1::T
    Δ::T
    lb::T
    ub::T
    fxn0::S
    fxn1::S
end

function _bisection(lb::Real, x::Real, ub::Real)
    if !(lb < x < ub)
        x = (ub+lb)*0.5
    end

    return x
end

function _update_brackets(lb::Real, fx::Real, x::Real, ub::Real)
    # cdf <= u
    fx ≤ 0 ? lb = x : ub = x

    return adjust_bracket((lb, ub))
end

function _newton_with_bisection(fxn0::Real, xn0::Real, F::Roots.Callable_Function, lb::Real, ub::Real, Δ::Real)
    lb, ub = _update_brackets(lb, fxn0, xn0, ub)

    xn1 = xn0 - Δ

    xn1 = _bisection(lb, xn1, ub)

    fxn1, Δ = F(xn1)

    return xn0, xn1, Δ, lb, ub, fxn0, fxn1
end

function Roots.init_state(::NewtonBracket, F::Roots.Callable_Function, x)
    # Starting point for Newton
    x0 = float(first(x[1]))
    lb = float(first(x[2]))
    ub = float(first(x[3]))

    fxn0, Δ = F(x0)
    xn0, xn1, Δ, lb, ub, fxn0, fxn1 = _newton_with_bisection(fxn0, x0, F, lb, ub, Δ)

    NewtonBracketState(xn0, xn1, Δ, lb, ub, fxn0, fxn1)
end

# bracketed Newton
function Roots.update_state(
    ::NewtonBracket,
    F,
    o::NewtonBracketState{T,S},
    options,
    l=NullTracks(),
) where {T,S}
    # try Newton, if out of bounds use bisection
    # Newton step
    xn0, xn1, Δ, lb, ub, fxn0, fxn1 = _newton_with_bisection(
        o.fxn1,
        o.xn1,
        F,
        o.lb,
        o.ub,
        o.Δ
    )

    @reset o.lb = lb
    @reset o.ub = ub

    @reset o.xn0 = o.xn1
    @reset o.xn1 = xn1
    @reset o.Δ = Δ
    @reset o.fxn0 = o.fxn1
    @reset o.fxn1 = fxn1

    return o, false
end