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
function continuousphase(signal, events::AbstractVector{Int}; centerfun=_center)
    signalcent = centerfun(signal, events)
    Hsignal = hilbert(signalcent)

    return angle.(Hsignal)
end

function _center(signal, events)
    mi, ma = avgextrema(signal, events)
    return signal .- Ref(mi + (ma - mi)/2)
end

function _centerd(signal, events)
    mi, ma = intervalextrema(signal, events)
    mi_μ, ma_μ = circmeand(mi), circmeand(ma)
    return signal .- Ref(mi_μ + (ma_μ - mi_μ)/2)
end

function relativephase(sigA, sigB, events; centerfun=_center)
    θA = continuousphase(sigA, events; centerfun)
    θB = continuousphase(sigB, events; centerfun)

    return unwrap(θA; range=2pi) - unwrap(θB; range=2pi)
end

function crpensemble(sigA, sigB, events; centerfun=_center)
    crp = relativephase(sigA, sigB, events; centerfun)

    ens_avg, ens_std = limitcycle(crp, events, mean=circmean, std=circstd)

    return rad2deg.(ens_avg), rad2deg.(ens_std)
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

function _linreg(x::AbstractVector{T}, y::AbstractVector{U}) where {T,U}
    V = promote_type(T, U)
    _linreg(convert(Vector{V}, collect(x)), convert(Vector{V}, collect(y)))
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
    β = _linreg(1:length(y),y)
    @. return y - (β[1] + β[2]*y)
end

"""
    detrend!(y)

Mutate `y` by removing its linear trend
"""
function detrend!(y)
    β = _linreg(1:length(y),y)
    @. y -= β[1] + β[2]*y
end

struct ForwardBackwardPad; end

"""
    centraldiff(x, order=1; dt=1, padding=NaN)

Take a numeric finite central difference along the first dimension of the array `x`, where
`dt` is the period between samples in `x`, and `padding` can be `nothing` (returning an
array which is 2 elements shorter than `x`), NaN, or `ForwardBackwardPad()` (which uses a
forward or backward finite difference for the first and last elements, respectively, of the
result).
"""
function centraldiff(x::AbstractVector{T}, order::Integer=1; dt=1, padding=ForwardBackwardPad()) where T
    if padding === nothing
        x′ = similar(x, T, length(x)-2)
        _centraldiff!(x′, x, Val(order), dt)
    else
        x′ = similar(x)
        _centraldiff!(x′, x, Val(order), dt, padding)
    end

    return x′
end

function centraldiff(
    x::Matrix{T}, order::Integer=1; dt=1, padding=ForwardBackwardPad()
) where T
    sz = size(x)
    dims=2
    if padding === nothing
        sz = ntuple(i -> i == 1 ? sz[1] - 2 : sz[i], length(sz))
        x′ = similar(x, sz)
        foreach(eachslice(x′; dims), eachslice(x; dims)) do vx′, vx
            _centraldiff!(vx′, vx, Val(order), dt)
        end
    else
        x′ = similar(x, sz)
        foreach(eachslice(x′; dims), eachslice(x; dims)) do vx′, vx
            _centraldiff!(vx′, vx, Val(order), dt, padding)
        end
    end

    return x′::Matrix{T}
end

@inline function _centraldiff!(x′::AbstractVector, x::AbstractVector, order, dt, padding)
    _centraldiff!(@view(x′[begin+1:end-1]), x, order, dt)
    x′[begin] = padding
    x′[end] = padding

    return nothing
end

@inline function _centraldiff!(x′::AbstractVector, x::AbstractVector, order::Val{1}, dt, ::ForwardBackwardPad)
    _centraldiff!(@view(x′[begin+1:end-1]), x, order, dt)
    x′[begin] = (x[begin+1] - x[begin]) / dt
    x′[end] = (x[end] - x[end-1]) / dt
end

@inline function _centraldiff!(x′::AbstractVector, x::AbstractVector, order::Val{2}, dt, ::ForwardBackwardPad)
    _centraldiff!(@view(x′[begin+1:end-1]), x, order, dt)
    x′[begin] = (x[begin+2] - 2*x[begin+1] + x[begin]) / dt^2
    x′[end] = (x[end] - 2*x[end-1] + x[end-2]) / dt^2

    return nothing
end

@inline function _centraldiff!(x′::AbstractVector, x::AbstractVector, ::Val{1}, dt)
    start = firstindex(x)
    last = lastindex(x) - 2
    @boundscheck (start, last) == (firstindex(x′), lastindex(x′)) || throw(BoundsError())

    @turbo for i in start:last
        x′[i] = (x[i+2] - x[i]) / 2dt
    end

    return nothing
end

@inline function _centraldiff!(x′::AbstractVector, x::AbstractVector, ::Val{2}, dt)
    start = firstindex(x)
    last = lastindex(x) - 2
    @boundscheck (start, last) == (firstindex(x′), lastindex(x′)) || throw(BoundsError())

    @turbo for i in start:last
        x′[i] = (x[i+2] - 2x[i+1] + x[i]) / dt^2
    end

    return nothing
end

