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

function Roots.init_state(::NewtonBracket, F::Roots.Callable_Function, x)
    # Starting point for Newton
    xn0 = float(first(x[1]))
    fxn0, Δ = F(xn0)
    xn1 = xn0 - Δ

    lb, ub = adjust_bracket(x[2:3])

    if !(lb < xn1 < ub)
        xn1 = (ub-lb)*0.5
    end

    fxn1, Δ = F(xn1)

    # cdf <= u
    if fxn1 <= 0
        lb = xn1
    else
        ub = xn1
    end

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
    t = o.xn1 - o.Δ

    # Bisection step
    if !(o.lb < t < o.ub)
        t = (o.ub-o.lb)*0.5
    end

    #Δ = cdf/pdf
    cdf, Δ = F(t)

    # cdf <= u
    if cdf <= 0
        @set! o.lb = t
    else
        @set! o.ub = t
    end
    lb, ub = adjust_bracket((o.lb, o.ub))
    @set! o.lb = lb
    @set! o.ub = ub

    @set! o.xn0 = o.xn1
    @set! o.xn1 = t
    @set! o.Δ = Δ
    @set! o.fxn0 = o.fxn1
    @set! o.fxn1 = cdf

    return o, false
end