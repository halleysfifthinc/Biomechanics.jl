"""
    timenormalize(data, events[, len])

Normalize the data (by column) to lengths of `len` bounded by the events.
"""
function timenormalize(data::AbstractArray, events::AbstractVector,len=100)
    cols = size(data,2)
    res = zeros((length(events)-1)*len, cols)

    fill_normdims!(res, data, events, len)

    return res
end

"""
    normalizeevents(nt1, nt2, events[, len])

Convert the events to normalized time between nt1 and nt2.
"""
function normalizeevents(nt1::T, nt2::T, events, len=100) where T
    nt1 < nt2 || throw(DomainError("nt1 must be less than nt2"))
    prod(nt1 .<= events .<= nt2) || throw(DomainError("events must be between normalizing events"))

    nt = normtime(nt1,nt2, len)
    return [ findfirst(x -> x >= ev, nt) for ev in events ]
end

"""
    normalizedeventtimes(nts, events[, len])

Convert the array of events to normalized time given the array of nts marking the
normalized event bounds.
"""
function normalizeevents(nts::AbstractVector{T}, events::AbstractVector{T}, len=100) where T
    nevents = similar(events)
    for i in 1:(length(nts)-1)
        rel = find(x -> nts[i] <= x <= nts[i+1], events)
        @views nevents[rel] .= normalizeevents(nts[i], nts[i+1], events[rel], len)
    end
    return nevents
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
    normtime(t1, t2[, len])

Returns a range [t1,t2) with a step of the difference of t2 and t1 divided by `len`.
`len` defaults to 100 points.
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
