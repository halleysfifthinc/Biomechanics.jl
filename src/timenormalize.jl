"""
    timenormalize(data, events::AbstractVector{Int}[, len=100])

Normalize the first dimension to lengths of `len` bounded by the events.
"""
function timenormalize(data::AbstractArray{T}, events::AbstractVector,len=100)
    dims = size(data)
    res = Array{T}(undef, (length(events)-1)*len, dims[2:end]...)

    fill_normdims!(res, data, events, len)

    return res
end

function fill_normdims!(res::AbstractArray,
                        data::AbstractArray,
                        events::Vector{Int},
                        len=100)
    (size(res,2) != size(data,2)) && throw(
        ArgumentError("`res` and `data` must have the same number of columns"))
    for i in 1:size(data,2)
        # Create interpolation object
        itp = interpolate(@view(data[:,i]), BSpline(Cubic(Line(Interpolations.OnGrid()))))

        fill_normstrides!(@view(res[:,i]), itp, events, len)
    end

    nothing
end

function fill_normstrides!(nstr::AbstractArray, str::AbstractArray, events::Vector,len=100)
    @assert mod(size(nstr,1),len) === 0
    @assert size(nstr,1) === (length(events)-1)*len

    # Create time normalized strides
    @inbounds for s in 1:length(events)-1
        nstr[(1:len).+(s-1)*len] = str(normtime(events[s], events[s+1], len))
    end

    nothing
end

"""
    timestoindices(t1::Real, t2::Real, times[, len=100])

Find the equivalent indices for `times` when the times `t1` and `t2` match the `begin` and
`end + 1` of an array of length `len`.

This is useful when an array has been timenormalized using `basetimes` as the normalizing events, but when the location of `times` in the normalized array is also needed.

# Examples
```jldoctest
julia> normedarray = timenormalize(collect(.01:.01:3), [132, 246])
100-element Vector{Float64}:
 1.32
 1.3314
 1.3428
 1.3541999999999996
[...]

julia> times = [1.56, 2.0]
2-element Vector{Float64}:
 1.56
 2.0

julia> idx = timestoindices(1.32, 2.46, times)
2-element Vector{Int64}:
 23
 61

julia> isapprox(normedarray[idx], times; atol=.02)
true
```
"""
function timestoindices(t1::T, t2::T, events, len=100) where T <: Real
    t1 < t2 || throw(DomainError("t1 must be less than t2"))
    all(t1 .<= events .<= t2) || throw(DomainError("events must be between normalizing events"))

    nt = normtime(t1, t2, len)
    return [ findfirst(≥(ev), nt) for ev in events ]
end

"""
    timestoindices(basetimes::AbstractVector{<:Real}, times[, len=100])

Find the equivalent indices for elements in `times` which are between `basetimes[i]` and
`basetimes[i+1]` which match the `begin` and `end + 1` of an array of length `len` for all successive pairs in `basetimes`.
"""
function timestoindices(basetimes::AbstractVector{T}, times::AbstractVector{T}, len=100) where T
    ntimes = similar(times)
    for i in 1:(length(nts)-1)
        rel = findall(x -> nts[i] ≤ x ≤ nts[i+1], times)
        @views ntimes[rel] .= timestoindices(nts[i], nts[i+1], times[rel], len)
    end
    return ntimes
end

"""
    normtime(t1, t2[, len=100])

Returns a range [t1,t2) with a length of `len`.
"""
function normtime(t1::T,t2::T,len=100) where T
    # st = (t2-t1)/len
    # t1:st:round(T, t2-st)
    range(t1, stop=t2, length=len+1)[1:len]
end

"""
    avgcycle(data, len=100)

Calculate the ensemble average and stdev of `data` assuming a normalized length of `len`
"""
function avgcycle(data::AbstractArray, len=100)
    mod(size(data,1),len) == 0 || throw(ArgumentError("length of data must be even"*
                                                      " multiple of len"))

    N, wid = size(data)
    normal = similar(data, (len, wid))
    vnormal = similar(data, (len, wid))

    @inbounds @views for dim in 1:wid
        for i in 1:len
            normal[i,dim] = mean(data[i:len:N, dim])
            vnormal[i,dim] = std(data[i:len:N, dim])
        end
    end

    return (normal, vnormal)
end
