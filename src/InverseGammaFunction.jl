module InverseGammaFunction

import Bessels
import SpecialFunctions: SpecialFunctions, loggamma, digamma
using Roots

export invgamma, invloggamma

# BigFloat values here are computed to prec = 256.
@inline _gamma(x::T) where {T <: AbstractFloat} = Bessels.gamma(x)
@inline _gamma(x::BigFloat) = SpecialFunctions.gamma(x)
@inline _get_xmin(::Type{T}) where {T <: AbstractFloat} = 1.4616321449683622
@inline _get_xmin(::Type{BigFloat}) = big"1.461632144968362341262659542325721328468196204006446351295988408598786440353795"
@inline _get_ymin(::Type{T}) where {T <: AbstractFloat} = 0.8856031944108887
@inline _get_ymin(::Type{BigFloat}) = big"0.8856031944108887002788159005825887332079515336699034488712001658751362274173978"

@inline _get_logymin(::Type{T}) where {T <: AbstractFloat} = -0.12148629053584961
@inline _get_logymin(::Type{BigFloat}) = big"-0.1214862905358496080955145571776915821513561731299990388637243731331352975757896"

# The Gamma function Γ: ℝ⁺ → ℝ⁺ has a minimum at
# ≈ 1.461632138429461, taking the value ≈ 0.8856031944108887
#   1.4616321449683623412626595423257
#
# The derivative of the gamma function is given by Γ'(x) = Γ(x) ψ⁽⁰⁾(x),
# where ψ⁽⁰⁾(x) is the digamma function.
"""
    invgamma(z; branch=:hi)

Return the inverse gamma function of `z`.

`z` must be a real number greater than 0.885603194410888 and less than `gamma(171)`.
If `branch` is `:hi`, then the returned value is in `(xmin, 171.0)` where `xmin=1.461632138429461`.
If `branch` is `:lo`, then the returned value is in `(0, xmin)`.
"""
function invgamma(y::T; branch::Symbol=:hi) where {T}
    # Minimum value of Γ(x)
    _min_x = _get_xmin(T)
    y > _get_ymin(T) ||
        throw(ArgumentError(lazy"Argument must be > 0.885603194410888. Got $y."))
    if branch === :hi
        xmax = T(171)
        xmin = _min_x
    elseif branch === :lo
        xmax = _min_x
        xmin = zero(T)
    else
        throw(ArgumentError(lazy"Keyword argument `branch` must be one of `:hi` or `:lo`"))
    end
    f = x -> _gamma(x) - y
    find_zero(f, (xmin, xmax), Roots.Chandrapatla())
end

# The low branch sometimes fails to converge
"""
    invloggamma(y; branch::Symbol=:hi)

Return the inverse loggamma function of `z`.

`z` must be a real number greater than -0.12148629053584961.
If `branch` is `:hi`, then the returned value is greater than `xmin=1.461632138429461`.
If `branch` is `:lo`, then the returned value is less than `xmin`.
"""
function invloggamma(y::T; branch::Symbol=:hi) where {T}
    y > _get_logymin(T) ||
        throw(ArgumentError(lazy"Argument must be > -0.12148629053584961. Got $y."))
    if branch === :hi
        x = y + 2
    elseif branch === :lo
        x = 10.0^(-y/2)
    else
        throw(ArgumentError(lazy"Keyword argument `branch` must be one of `:hi` or `:lo`"))
    end
    i = 1
    while true
        incr = (loggamma(x) - y) / digamma(x)
        abs(incr) <= 3 * eps(x) && return x
        x = x - incr
        i += 1
        i > 1000 && break
    end
    @show (loggamma(x) - y) / digamma(x)
    @show i, x
    throw(ArgumentError(lazy"invloggamma($y) failed to converge"))
end

end # module InverseGammaFunction
