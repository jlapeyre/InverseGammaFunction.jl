module InverseGammaFunction

using SpecialFunctions: gamma
using Roots

export invgamma

"""
    invgamma(z; branch=:hi)

Return the inverse gamma function of `z`.

`z` must be a real number greater than 0.885603194410888 and less than `gamma(171)`.
If `branch` is `:hi`, then the returned value is in `(xmin, 171.0)` where `xmin=1.461632138429461`.
If `branch` is `:lo`, then the returned value is in `(0, xmin)`.
"""
function invgamma(z; branch=:hi)
    z >= 0.8856031944108887 ||
        throw(ArgumentError(lazy"Argument must be > 0.885603194410888. Got $z."))
    if branch == :hi
        xmax = 171.0
        xmin = 1.461632138429461
    elseif branch == :lo
        xmax = 1.461632138429461
        xmin = 0.0
    else
        throw(ArgumentError(lazy"Keyword argument `branch` must be one of `:hi` or `:lo`"))
    end
    f = x -> gamma(x) - z
    find_zero(f, (xmin, xmax), Roots.Chandrapatla())
end

end # module InverseGammaFunction
