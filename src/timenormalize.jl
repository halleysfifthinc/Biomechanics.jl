"""
    timenormalize(data, events::Vector{Int}, [len=100])
    timenormalize(data, intervals::Vector{<:AbstractRange{Int}}, [len=100])

Normalize the first dimension to lengths of `len` bounded by the events.
"""
function timenormalize(data::AbstractArray{T}, events, len::Int=100) where T
    dims = size(data)
    oobevents = events .∈ Ref(axes(data,1))
    any(~, oobevents) &&
        throw(error("events $(events[.~oobevents]) are not valid indices of `data`"))
    res = Array{T}(undef, (length(events)-1)*len, dims[2:end]...)

    fill_normdims!(res, data, intervals(events), len)

    return res
end

⊂ = issubset
function timenormalize(data::AbstractArray{T}, intvls::Vector{<:AbstractRange{Int}}, len::Int=100) where T
    dims = size(data)
    oobevents = intvls .⊂ Ref(axes(data,1))
    any(~, oobevents) &&
        throw(error("intervals $(intvls[.~oobevents]) are not valid indices of `data`"))
    res = Array{T}(undef, length(intvls)*len, dims[2:end]...)

    fill_normdims!(res, data, intvls, len)

    return res
end

function fill_normdims!(res, data, intvls, len::Int=100)
    (size(res,2) != size(data,2)) && throw(
        ArgumentError("`res` and `data` must have the same number of columns"))
    for i in 1:size(data,2)
        # Create interpolation object
        itp = interpolate(@view(data[:,i]), BSpline(Cubic(Line(Interpolations.OnGrid()))))

        fill_normstrides!(@view(res[:,i]), itp, intvls, len)
    end

    nothing
end

function fill_normstrides!(nstr, str, intvls::Vector{<:AbstractRange{Int}}, len::Int=100)
    @assert mod(size(nstr,1),len) === 0
    @assert size(nstr,1) === length(intvls)*len

    # Create time normalized strides
    @inbounds for s in eachindex(intvls)
        nstr[(1:len).+(s-1)*len] = str(normtime(intvls[s], len))
    end

    nothing
end

function time_denormalize(data, events, len::Int=100; fs)
    @views normed_time = mapreduce((a,b) -> normtime(a,b,len), vcat,
        events[begin:end-1], events[begin+1:end])
    interp = linear_interpolation(normed_time, data)
    denormed_time = first(events):inv(fs):last(normed_time)

    return interp.(denormed_time)
end

"""
    timestoindices(t1::Real, t2::Real, times[, len=100])

Find the equivalent indices for `times` in a range `[t1,t2)` of length `len`.

Equivalent indices are calculated as the index of the first element (in `[t1,t2)`) with a
value greater than or after the corresponding value in `times`.

This is useful when an array has been timenormalized using one set of events for
normalizing, but another set of events needs to be used in indexing into the timenormalized
array.

# Examples
```jldoctest
julia> normedarray = collect(normtime(1.32, 2.46, 100))
100-element Vector{Float64}:
 1.32
 1.3314
 1.3428
 1.3542
[...]

julia> times = [1.56, 2.0];

julia> idx = timestoindices(1.32, 2.46, times)
2-element Vector{Int64}:
 23
 61

julia> normedarray[idx]
2-element Vector{Float64}:
 1.5708
 2.004

```
"""
function timestoindices(t1::T, t2::T, events, len::Int=100) where T <: Real
    t1 < t2 || throw(DomainError("t1 must be less than t2"))
    isempty(events) && return Int[]
    all(t1 .<= events .<= t2) || throw(DomainError("events must be between normalizing events"))

    nt = normtime(t1, t2, len)
    return Int[ searchsortedfirst(nt, ev) for ev in events ]
end

"""
    timestoindices(basetimes::AbstractVector{<:Real}, times[, len=100])

Find the equivalent indices for elements in `times` which are in a range `[t1,t2)` of length
`len`, where `t1` and `t2` are `basetimes[i]` and `basetimes[i+1]`, for all successive pairs
in `basetimes`.

Equivalent indices are calculated as the index of the first element (in `[t1,t2)`) with a
value greater than or after the corresponding value in `times`.

# Examples
```jldoctest
julia> normedarray = collect(normtime(1.32, 2.46, 100))
100-element Vector{Float64}:
 1.32
 1.3314
 1.3428
 1.3542
[...]

julia> times = [1.56, 2.0];

julia> idx = timestoindices([1.32, 2.46], times)
2-element Vector{Int64}:
 23
 61

julia> normedarray[idx]
2-element Vector{Float64}:
 1.5708
 2.004

```
"""
function timestoindices(
    basetimes::AbstractVector{T}, times::AbstractVector{T}, len::Int=100
) where T
    ntimes = similar(times, Int)
    perm = sortperm(times)
    stimes = permute!(copy(times), perm)
    for i in 1:(length(basetimes)-1)
        f = searchsortedfirst(stimes, basetimes[i])
        l = searchsortedlast(stimes, basetimes[i+1])
        @views ntimes[f:l] .= timestoindices(basetimes[i], basetimes[i+1], stimes[f:l], len)
    end

    return invpermute!(ntimes, perm)
end

export timeintervals_toindices
"""
    timeintervals_toindices(t1::Real, t2::Real, ranges[, len=100])

Find the equivalent ranges of indices in a range `[t1,t2)` of length `len` for each range in
`ranges` which overlap `[t1,t2)`. Only the indices for the overlapping portion will be
returned for a range which only partially overlaps `[t1,t2)`.

(Equivalent to `timestoindices`, but for ranges instead of scalars.)

See also: [`timestoindices`](@ref)

# Examples
```jldoctest
julia> normedarray = collect(normtime(1.32, 2.46, 100))
100-element Vector{Float64}:
 1.32
 1.3314
 1.3428
 1.3542
[...]

julia> times = [1.56, 2.0];

julia> idx = timeintervals_toindices(1.32, 2.46, 1.56:2.0)
2-element Vector{Int64}:
 23
 61

julia> normedarray[idx]
2-element Vector{Float64}:
 1.5708
 2.004

```
"""
function timeintervals_toindices(
    t1::T, t2::T, ranges::Vector{<:AbstractRange}, len::Int=100
) where T
    t1 < t2 || throw(DomainError("t1 must be less than t2"))
    ranges = sort(ranges)

    nt = normtime(t2, t2, len)
    return [
        searchsortedfirst(nt, first(ev)):searchsortedlast(nt, last(ev))
            for ev in ranges
        ]
end

function timeintervals_toindices(
    t1::T, t2::T, rg::AbstractRange, len::Int=100
) where T
    return timeintervals_toindices(t1, t2, [rg], len)
end

function timeintervals_toindices(basetimes::Vector{T}, events::Vector{<:AbstractRange}, len::Int=100) where T
    timeintervals_toindices(intervals(basetimes; endincluded=true, step=nothing), events, len)
end

function timeintervals_toindices(baseranges::Vector{<:AbstractRange}, events::Vector{<:AbstractRange}, len::Int=100)
    events = sort(events)

    ntimes = UnitRange[]
    for baserange in baseranges
        overlapped = findall(ev -> !isdisjoint(baserange, ev), events)
        append!(ntimes, timeintervals_toindices(first(baserange), last(baserange), events[overlapped], len))
    end

    return ntimes
end

"""
    normtime(t1, t2[, len=100])

Returns a range [t1,t2) with a length of `len`.
"""
function normtime(t1, t2, len::Int=100)
    # st = (t2-t1)/len
    # t1:st:round(T, t2-st)
    range(t1, stop=t2, length=len+1)[1:len]
end

function normtime(rg, len::Int=100)
    range(rg.start, stop=rg.stop, length=len+1)[1:len]
end

"""
    limitcycle(data, [events::AbstractVector{Int}, len=100])

Calculate the ensemble average and stdev of `data` with a normalized length of `len`
"""
function limitcycle(data, len::Int=100; mean=mean, std=std)
    N = size(data, 1)
    mod(N, len) == 0 ||
        throw(ArgumentError("length of data must be even multiple of len"))

    wid = size(data, 2)
    ensemble_avg = similar(data, (len, wid))
    ensemble_std = similar(data, (len, wid))

    @inbounds @views for dim in 1:wid
        for i in 1:len
            ensemble_avg[i,dim] = mean(data[i:len:N, dim])
            ensemble_std[i,dim] = std(data[i:len:N, dim])
        end
    end

    return (vec(ensemble_avg), vec(ensemble_std))
end

function limitcycle(data, events, len::Int=100; mean=mean, std=std)
    normed = timenormalize(data, events, len)
    limitcycle(normed, len; mean, std)
end
