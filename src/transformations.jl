"""
    continuousphase(signal, events) -> Vector

Calculate the continuous phase of `signal` by finding the angle of the Hilbert transform of the
amplitude centered signal, `signal`. The amplitude centering uses the average extrema between
pairs of `events`.

This method of continuous phase calculation is as recommended by Lamb and Stöckl (2014) [1]

[1] P. F. Lamb and M. Stöckl, “On the use of continuous relative phase: Review of current
    approaches and outline for a new standard,” Clinical Biomechanics, vol. 29, no. 5,
    pp. 484–493, May 2014, [doi](https://doi.org/10.1016%2Fj.clinbiomech.2014.03.008).
"""
function continuousphase(signal, events::AbstractVector{Int})
    mi, ma = avgextrema(signal, events)
    signalcent = signal .- Ref(mi - (ma - mi)/2)
    Hsignal = hilbert(signalcent)

    return angle.(Hsignal)
end

"""
    demean(x)

Subtract the mean of `x` from `x`
"""
function demean(x)
    return x .- mean(x)
end

"""
    demean!(x)

Mutate `x` by subtracting its mean
"""
function demean!(x)
    x .-= mean(x)
end

# From @dmbates in https://discourse.julialang.org/t/efficient-way-of-doing-linear-regression/31232/28
function _linreg(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:AbstractFloat}
    (N = length(x)) == length(y) || throw(DimensionMismatch())
    ldiv!(cholesky!(Symmetric([T(N) sum(x); zero(T) sum(abs2, x)], :U)), [sum(y), dot(x, y)])
end

"""
    detrend(y)

Remove the linear trend from `y`
"""
function detrend(y)
    β = linreg(1:length(y),y)
    @. return y - (β[1] + β[2]*y)
end

"""
    detrend!(y)

Mutate `y` by removing its linear trend
"""
function detrend!(y)
    a, b = linreg(1:length(y),y)
    @. y -= β[1] + β[2]*y
end

struct ForwardBackwardPad; end

"""
    centraldiff(x, order; dt=1, padding=NaN)

Take a numeric finite central difference of the array `x`, where `dt` is the period between
samples in `x`, and `padding` can be `nothing` (returning an array which is 2 elements
shorter than `x`), NaN, or `ForwardBackwardPad()` (which uses a forward or backward finite
difference for the first and last elements, respectively, of the result).
"""
function centraldiff(x::AbstractVector{T}, order::Integer; dt=1, padding=NaN) where T
    if padding === nothing
        x′ = similar(x, T, length(x)-2)
        _centraldiff!(x′, x, Val(order), dt)
    else
        x′ = similar(x)
        _centraldiff!(x′, x, Val(order), dt, padding)
    end

    return x′
end

@inline function _centraldiff!(x′, x, order::Val{1}, dt)
    start = firstindex(x)
    last = lastindex(x) - 2
    @boundscheck (start, last) == (firstindex(x′), lastindex(x′)) || throw(BoundsError())

    @inbounds @simd ivdep for i in start:last-2
        x′[i] = (x[i+2] - x[i]) / 2dt
    end

    return nothing
end

@inline function _centraldiff!(x′, x, order::Val{2}, dt)
    start = firstindex(x)
    last = lastindex(x) - 2
    @boundscheck (start, last) == (firstindex(x′), lastindex(x′)) || throw(BoundsError())

    @inbounds @simd ivdep for i in start:last-2
        x′[i] = (x[i+2] - 2x[i+1] + x[i]) / dt^2
    end

    return nothing
end

@inline function _centraldiff!(x′, x, order, dt, padding)
    _centraldiff!(@view(x′[begin+1:end-1]), x, order, dt)
    x′[begin] = padding
    x′[end] = padding

    return nothing
end

@inline function _centraldiff!(x′, x, order::Val{1}, dt, padding::ForwardBackwardPad)
    _centraldiff!(@view(x′[begin+1:end-1]), x, order, dt)
    x′[begin] = (x[begin+1] - x[begin]) / dt
    x′[end] = (x[end] - x[end-1]) / dt
end

@inline function _centraldiff!(x′, x, order::Val{2}, dt, padding::ForwardBackwardPad)
    _centraldiff!(@view(x′[begin+1:end-1]), x, order, dt)
    x′[begin] = (x[begin+2] - 2*x[begin+1] + x[begin]) / dt^2
    x′[end] = (x[end] - 2*x[end-1] + x[end-2]) / dt^2

    return nothing
end

