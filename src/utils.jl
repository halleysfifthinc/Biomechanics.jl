"""
    intervals(a, [b]; endclosed=false, step=one(eltype(a)))

Return a vector of ranges from `a[i]:step:a[i+1]` or `a[i]:step:b[j]`, if `b` is given.

All ranges returned will have a non-zero length and will not overlap (except by one `step` if endclosed is `true`). Values in `a` and/or `b` will be skipped when required to enforce this behavior.

If `endclosed == true`, ranges will include the `a[i+1]` or `b[j]` value, otherwise, ranges will end one `step` before that value.

"""
function intervals(a::AbstractVector{T}; endclosed=false, step=one(T)) where T
    length(a) > 1 || throw(ArgumentError("a is not long enough to create any intervals"))
    ranges = Vector{typeof(first(a):step:a[2])}()
    for i in 1:(length(a)-1)
        if endclosed
            rg = a[i]:step:a[i+1]
        else
            rg = (a[i]:step:a[i+1])[1:end-1]
        end
        push!(ranges, rg)
    end

    return ranges
end

function intervals(
    a::AbstractVector{T}, b::AbstractVector{T}; endclosed=false, step=one(T)
) where T
    ranges = Vector{typeof(first(a):step:first(b))}()

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

function _rotating_diff(vectors::Vararg{<:AbstractVector{T},N};
    i=argmin(first.(vectors))
) where {T,N}
    diffed = ntuple(x -> similar(Vector{T}, (0,)), N)
    iterators::Vector{Union{Nothing,Tuple{T,Any}}} = collect(iterate.(vectors))

    nextvector = i
    currvector = nextvector
    nextvector = mod1(nextvector+1, N)

    while !isnothing(iterators[nextvector])
        cval, cstate = iterators[currvector]
        nval, nstate = iterators[nextvector]
        iterators[currvector] = iterate(vectors[currvector], cstate)

        push!(diffed[currvector], nval - cval)
        currvector = nextvector
        nextvector = mod1(nextvector+1, N)
    end

    return diffed
end

function interleave(vectors::Vararg{<:AbstractVector{T},N}) where {T,N}
    minlen = minimum(length.(vectors))
    arr = Vector{T}(undef, N*minlen)
    iterators::Vector{Union{Nothing,Tuple{T,Any}}} = collect(iterate.(vectors))

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

