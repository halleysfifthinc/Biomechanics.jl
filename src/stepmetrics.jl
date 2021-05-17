"""
    SpatiotemporalMetrics

Functions for calculating spatiotemporal gait characteristics.

FS = Footstrike
FO = Foot liftoff

Stride time = time between FS of the same foot
Step time = time since previous opposite FS (e.g. left step time is the time between the current left step and the most recent right step)

Swing time = time from FO to FS of the same leg
Stance time = time from FS to FO of the same leg

Single support = time when only one foot is in contact with the ground; equivalent to the
sum of left and right swing when walking, and the sum of the left and right stance when
running

Double support = time when both feet are in contact with the ground; equivalent to the stride time minus the sum of the left and right swing when walking, and zero when running

Walking:

|RFS|     |LFO|     |LFS|     |RFO|     |RFS|     |LFO|     |LFS|
-----------------------------------------------------------------
  |         Right stance        | R swing |
  ---------------------------------------------------
            | L swing |            Left stance      |
  ---------------------------------------------------
  |  L DS   |         |  R DS   |         |
  -----------------------------------------
  |  Left step time   |  Right step time  |
  -----------------------------------------
  |             Stride time               |
  -----------------------------------------

Running:

|RFS|      |RFO|     |LFS|      |LFO|     |RFS|     |RFO|     |LFS|
-------------------------------------------------------------------
  | R stance |         Right swing          |  Float  |
  ---------------------------------------------------------------
             |  Float  | L stance |           Left swing        |
  ---------------------------------------------------------------
  |  Left step time    |  Right step time   |
  -------------------------------------------
  |             Stride time                 |
  -------------------------------------------
"""
module SpatiotemporalMetrics

using ..Biomechanics, LinearAlgebra
using ..Biomechanics: _rotating_diff

export beginwithevent, stridetimes, steptimes, swingstance, swing, stance, singlesupport,
    doublesupport, steplength, stepwidth

function beginwithevent(fs, events...)
    ffs = something(findlast(<(maximum(first, events)), fs), 1)
    pred = >(fs[ffs])
    fi = findfirst.(pred, events)

    return (@view(fs[ffs:end]), view.(events, range.(fi, lastindex.(events); step=1))...)
end

"""
    steplength(;lfs, rfs, lftpos, rftpos, [AP=1, VT=3, fs=1, ref=:fixed]) -> (;lsteps, rsteps)

Calculate left and right step lengths.

If the eltype's of `lfs` or `rfs` are Integers, they are assumed to be already in units of samples.

# Keyword Arguments:
- `AP=1`: The column for the anteroposterior axis
- `VT=1`: The column for the vertical axis
- `fs=1`: The sampling frequency
"""
function steplength(;lfs, lfo, rfs, rfo, lftpos, rftpos, AP=1, VT=3, fs=1)
    size(lftpos, 1) == size(rftpos, 1) ||
        throw(ArgumentError("left and right foot position data have unequal lengths"))

    if eltype(lfs) <: AbstractFloat
        lfs = round.(Int, lfs .* fs)
    end
    if eltype(lfo) <: AbstractFloat
        lfo = round.(Int, lfo .* fs)
    end
    if eltype(rfs) <: AbstractFloat
        rfs = round.(Int, rfs .* fs)
    end
    if eltype(rfo) <: AbstractFloat
        rfo = round.(Int, rfo .* fs)
    end

    max(last(lfs), last(rfs)) ≤ size(lftpos, 1) ||
        throw(ArgumentError("foot strike after the end of foot position data"))

    steplength(lfs, lfo, rfs, rfo; lftpos, rftpos, AP, VT)
end

function steplength(
    lfs::AbstractVector{Int},
    lfo::AbstractVector{Int},
    rfs::AbstractVector{Int},
    rfo::AbstractVector{Int};
    lftpos, rftpos, AP=1, VT=3
)
    size(lftpos, 1) == size(rftpos, 1) ||
        throw(ArgumentError("left and right foot position data have unequal lengths"))
    max(last(lfs), last(rfs)) ≤ size(lftpos, 1) ||
        throw(ArgumentError("foot strike after the end of foot position data"))
    any(iszero, doublesupport(;lfs, lfo, rfs, rfo)) &&
        throw(ArgumentError("running step (double support = 0) encountered. this function only support data from walking"))

    if VT === nothing
        slc = [AP]
    else
        slc = sort([AP,VT])
    end

    rsteps = vec(mapslices(norm, rftpos[rfs, slc] - lftpos[rfs, slc]; dims=2))
    lsteps = vec(mapslices(norm, lftpos[lfs, slc] - rftpos[lfs, slc]; dims=2))

    return (;lsteps, rsteps)
end

function stepwidth(;lfs, rfs, lftpos, rftpos, ML=2, fs=1)
    size(lftpos, 1) == size(rftpos, 1) ||
        throw(ArgumentError("left and right foot position data have unequal lengths"))

    lfs = round.(Int, lfs .* fs)
    rfs = round.(Int, rfs .* fs)

    max(last(lfs), last(rfs)) ≤ size(lftpos, 1) ||
        throw(ArgumentError("foot strike after the end of foot position data"))

    stepwidth(lfs, rfs; lftpos, rftpos, AP, ML)
end

function stepwidth(lfs::AbstractVector{Int}, rfs::AbstractVector{Int};
    lftpos, rftpos, ML=2
)
    size(lftpos, 1) == size(rftpos, 1) ||
        throw(ArgumentError("left and right foot position data have unequal lengths"))
    max(last(lfs), last(rfs)) ≤ size(lftpos, 1) ||
        throw(ArgumentError("foot strike after the end of foot position data"))

    if first(lfs) > first(rfs)
        lsteps, rsteps = _rotating_diff(rftpos[rfs, ML], lftpos[lfs, ML]; i=1)
    else
        rsteps, lsteps = _rotating_diff(lftpos[lfs, ML], rftpos[rfs, ML]; i=1)
    end

    return (;lsteps, rsteps)
end

stridetimes(fs) = diff(fs)

"""
    steptimes(;lfs, rfs) -> (;lsteps, rsteps)

Calculate the left and right step times, where a left step is the left foot-strike
following a right foot-strike, and vice versa.
"""
function steptimes(;lfs, rfs)
    f = findfirst(>(first(lfs)), rfs)
    if first(lfs) > first(rfs)
        lsteps, rsteps = _rotating_diff(rfs, lfs)
    else
        rsteps, lsteps = _rotating_diff(lfs, rfs)
    end

    return (;lsteps, rsteps)
end

"""
    swingstance(fs, fo; normalize=true) -> (;swing, stance)

Calculate the swing and stance times beginning with the first fs.

If `normalize` is true, swing and stance times will be normalized to the stride time given by `fs`. This will strip any `fo` preceding the first `fs` and any `fo` following the last `fs`. (This would remove an initial swing and/or final stance.)
"""
function swingstance(fs, fo; normalize=true)
    if normalize
        fs, fo = beginwithevent(fs, fo)
    end

    if first(fo) < first(fs)
        sw, st = _rotating_diff(fo, fs)
    else
        st, sw = _rotating_diff(fs, fo)
    end

    if normalize
        l = min(length(st), length(sw))
        resize!(st, l)
        resize!(sw, l)

        strides = stridetimes(fs)
        sw ./= strides
        st ./= strides
    end

    return (;swing=sw, stance=st)
end

function swing(fs, fo; normalize=true)
    sw, _ = swingstance(fs, fo; normalize)
    return sw
end

function stance(fs, fo; normalize=true)
    _, st = swingstance(fs, fo; normalize)
    return st
end

"""
    singlesupport(lfs, lfo, rfs, rfo)

Calculate the single support times based on right strides and beginning with the first rfs.

See also: [`swingstance`](@ref), [`swing`](@ref), [`stance`](@ref)
"""
function singlesupport(;lfs, lfo, rfs, rfo, normalize=true)
    rfs, lfo, lfs, rfo = beginwithevent(rfs, lfo, lfs, rfo)

    lsteps, rsteps = steptimes(;rfs,lfs)
    lswings, lstances = swingstance(lfs, lfo; normalize=false)
    rswings, rstances = swingstance(rfs, rfo; normalize=false)
    l = minimum(lastindex, (lsteps, rsteps, lswings, lstances, rswings, rstances))

    ss = @. @views ifelse(rstances[1:l] > lsteps[1:l], rswings[1:l], rstances[1:l]) +
        ifelse(lstances[1:l] > rsteps[1:l], lswings[1:l], lstances[1:l])

    any(≤(0), ss) && throw(error("Single support cannot be ≤ 0. Check your gait events"))

    if normalize
        strides = stridetimes(rfs)
        ss ./= strides[1:l]
    end

    return ss
end

"""
    doublesupport(lfs, lfo, rfs, rfo)

Calculate the double support times based on right strides and beginning with the first rfs.

See also: [`swingstance`](@ref), [`swing`](@ref), [`stance`](@ref)
"""
function doublesupport(;lfs, lfo, rfs, rfo, normalize=true)
    rfs, lfo, lfs, rfo = beginwithevent(rfs, lfo, lfs, rfo)
    lsteps, rsteps = steptimes(;rfs,lfs)
    lswings, lstances = swingstance(lfs, lfo; normalize=false)
    rswings, rstances = swingstance(rfs, rfo; normalize=false)
    l = minimum(lastindex, (lsteps, rsteps, lswings, lstances, rswings, rstances))

    ds = @. @views ifelse(rstances[1:l] < lsteps[1:l], (0,), lsteps[1:l] - lswings[1:l]) +
        ifelse(lstances[1:l] < rsteps[1:l], (0,), rsteps[1:l] - rswings[1:l])

    if normalize
        strides = stridetimes(rfs)
        ds = ds ./ strides[1:l]
    end

    return ds
end

function swingintervals(fs::AbstractVector{T}, fo::AbstractVector{T}; step=one(T)) where T
    intervals(fo, fs; step)
end

function stanceintervals(fs::AbstractVector{T}, fo::AbstractVector{T}; step=one(T)) where T
    intervals(fs, fo; step)
end

end
