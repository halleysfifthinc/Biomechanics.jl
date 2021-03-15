stridetime(fs) = diff(fs)

"""
    steptimes(lfs, rfs) -> (leftsteps, rightsteps)

Calculate the left and right step times, where a left step is the left foot-strike
following a right foot-strike, and vice versa.
"""
steptimes(lfs, rfs) = _rotating_diff(rfs, lfs)

"""
    swingstance(fs, fo) -> (swing, stance)

Calculate the swing and stance times beginning with the first fs.
"""
function swingstance(fs, fo)
    f = findfirst(>(first(fs)), fo)
    sw, st = _rotating_diff(@view(fo[f:end]), fs)
end

function swing(fs, fo)
    sw, _ = swingstance(fs, fo)
    return sw
end

function stance(fs, fo)
    _, st = swingstance(fs, fo)
    return st
end

function intervals(a::AbstractVector{T}, b::AbstractVector{T}; step=one(T)) where T
    ranges = Vector{typeof(first(a):first(b))}()
    sizehint!(ranges, length(a))

    firstbi = searchsortedfirst(b, first(a))

    bi = firstbi
    ai = firstindex(a)
    while bi < lastindex(b) && ai < lastindex(a)
        bi = searchsortedfirst(b, a[ai], bi, lastindex(b), Base.Order.Forward)
        ai = searchsortedlast(a, b[bi], ai, lastindex(a), Base.Order.Forward)
        push!(ranges, a[ai]:step:b[bi])

        bi += 1
        ai += 1
    end

    return ranges
end

function swingintervals(fs::AbstractVector{T}, fo::AbstractVector{T}; step=one(T)) where T
    intervals(fo, fs; step)
end

function stanceintervals(fs::AbstractVector{T}, fo::AbstractVector{T}; step=one(T)) where T
    intervals(fs, fo; step)
end

"""
    singlesupport(lfs, lfo, rfs, rfo)

Calculate the single support times based on right strides and beginning with the first rfs.

See also: [`swingstance`](@ref), [`swing`](@ref), [`stance`](@ref)
"""
function singlesupport(lfs, lfo, rfs, rfo)
    flfo = findf(>(f(rfs)), lfo)
    frfo = findf(>(f(rfs)), rfo)
    flfs = findf(>(f(rfs)), lfs)
    @views _, leftswing, _, rightswing = _rotating_diff(rfs, lfo[flfo:end], lfs[flfs:end], rfo[frfo])

    l = min(lastindex(leftswing), lastindex(rightswing))
    @views ss = leftswing[1:l] + rightswing[1:l]

    return ss
end

"""
    doublesupport(lfs, lfo, rfs, rfo)

Calculate the double support times based on right strides and beginning with the first rfs.

See also: [`swingstance`](@ref), [`swing`](@ref), [`stance`](@ref)
"""
function singlesupport(lfs, lfo, rfs, rfo)
    ss = singlesupport(lfs, lfo, rfs, rfo)
    strides = stridetimes(rfs)

    l = min(lastindex(strides), lastindex(ss))
    @views ds = strides[1:l] - ss[1:l]

    return ds
end

"""
    stepsingledoublesupport(lfs, lfo, rfs, rfo) -> (rds, lss, lds, rss)

Calculate the right double-support, left single-support (aka swing), left double-support,
and right single-support (aka swing).
"""
function singledoublesupport(lfs, lfo, rfs, rfo)
    rds, lss, lds, rss = _rotating_diff(rfs, lfo, lfs, rfo)
end

function _rotating_diff(vectors::Vararg{Vector{T},N}) where {T,N}
    diffed = ntuple(x -> similar(Vector{T}, (0,)), N)
    iterators::Vector{Union{Nothing,Tuple{T,Int}}} = collect(iterate.(vectors))

    nextvector = argmin(first.(vectors))
    currvector = nextvector
    nextvector = mod1(nextvector+1, N)

    while !isnothing(iterators[nextvector])
        cval, cstate = iterators[currvector]
        nval, nstate = iterators[nextvector]
        iterators[currvector] = iterate(vectors[currvector], cstate)

        push!(diffed[currvector], round(nval - cval; digits=2))
        currvector = nextvector
        nextvector = mod1(nextvector+1, N)
    end

    return diffed
end

function interleave(vectors::Vararg{<:AbstractVector{T},N}) where {T,N}
    minlen = minimum(length.(vectors))
    arr = Vector{T}(undef, N*minlen)
    iterators::Vector{Union{Nothing,Tuple{T,Int}}} = collect(iterate.(vectors))

    currvector = 1
    i = 0

    while !isnothing(iterators[currvector])
        cval, cstate = iterators[currvector]
        iterators[currvector] = iterate(vectors[currvector], cstate)

        arr[i += 1] = cval
        currvector = mod1(currvector+1, N)
    end

    return arr
end

function circshift1(A,shift)
    l = length(A)
    if abs(shift) == l
        return A
    else
        return (A[mod1(1+shift,l):end]...,A[1:mod1(shift,l)]...)
    end
end

function circshift2(A,shift)
    l = length(A)
    if abs(shift) == l
        return A
    else
        return ntuple(i -> A[mod1(i+shift, l)], l)
    end
end

