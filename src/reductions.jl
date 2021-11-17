"""
    circmean(x)

Compute the circular mean of `x` in radians.

[1] N. I. Fisher, Statistical Analysis of Circular Data. Cambridge University Press, 1993.
"""
function circmean(x)
    s = mean(sin, x)
    c = mean(cos, x)

    atan(s, c)
end

"""
    circmeand(x)

Compute the circular mean of `x` in degrees.

[1] N. I. Fisher, Statistical Analysis of Circular Data. Cambridge University Press, 1993.
"""
function circmeand(x)
    s = mean(sind, x)
    c = mean(cosd, x)

    atand(s, c)
end

"""
    circstd(x)

Compute the circular standard deviation of `x` in radians.

[1] N. I. Fisher, Statistical Analysis of Circular Data. Cambridge University Press, 1993.
"""
function circstd(x)
    s = mean(sin, x)
    c = mean(cos, x)

    sqrt(-2*log(hypot(c, s)))
end

"""
    circstdd(x)

Compute the circular standard deviation of `x` in degrees.

[1] N. I. Fisher, Statistical Analysis of Circular Data. Cambridge University Press, 1993.
"""
function circstdd(x)
    s = mean(sind, x)
    c = mean(cosd, x)

    rad2deg(sqrt(-2*log(hypot(c, s))))
end

"""
    mean_std_range(x, events) -> Tuple(avg, std)

Find the average range and range variability of `x` for the set of all intervals given by
`events`.
"""
function mean_std_range(x, events::AbstractVector{Int})
    mi, mx = intervalextrema(x, events)
    roms = mx - mi

    return mean(roms), std(roms)
end

"""
    avgextrema(x, events) -> Tuple(min, max)

Find the average minima and maxima for the set of all intervals given by `events`.
"""
function avgextrema(x, events::AbstractVector{Int})
    mi, mx = intervalextrema(x, events)

    return mean(mi), mean(mx)
end

function intervalextrema(x::AbstractVector{T}, events) where T
    rgs = intervals(events; endincluded=false)
    mi = Vector{T}(undef, length(rgs))
    mx = Vector{T}(undef, length(rgs))

    intervalextrema!(mi, mx, x, rgs)
end

function intervalextrema!(mi, mx, x, rgs)
    @inbounds for (i, rg) in enumerate(rgs)
        mi[i], mx[i] = extrema(view(x, rg))
    end

    return mi, mx
end

