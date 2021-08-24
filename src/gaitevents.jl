function toindices(times, fs)
    round.(Int, times.*fs) .+ 1
end

function totimes(indices, fs)
    (indices .- 1)./fs
end

"""
    findduplicates(itr)

"""
function findduplicates(itr)
    dups = Dict{eltype(itr), Vector{Int}}()

    for (i, v) in enumerate(itr)
        push!(get!(dups, v, Int[]), i)
    end
    filter!(((k,v),) -> length(v) !== 1, dups)

    return dups
end

function matchevents(pred::AbstractVector, act::AbstractVector; Tthresh=0.5*median(diff(act)))
    mdiff_pred = median(diff(pred))
    mdiff_act = median(diff(act))
    err(a,b) = abs(a - b)/b
    if err(mdiff_pred, mdiff_act) > 10 || err(mdiff_act, mdiff_pred) > 10
        throw(error("""detected a large difference in event frequency. Are predicted events
            in the same units as actual events (eg frames vs sec)?"""))
    end
    missed = max(0, length(act) - length(pred))

    _pred = copy(pred)
    tree = KDTree(reshape(act, 1, length(act)), Chebyshev())
    idxs, dists = nn(tree, reshape(pred, 1, length(pred)))

    dups = findduplicates(idxs)
    delidxs = Int[]
    if !isempty(dups)
        foreach(dups) do (k, v)
            append!(delidxs, v[setdiff(eachindex(v), argmin(dists[v]))])
        end
        sort!(delidxs)
        deleteat!(idxs, delidxs)
        deleteat!(dists, delidxs)
        deleteat!(_pred, delidxs)
    end

    dists .= _pred - act[idxs]

    return setdiff(eachindex(pred), delidxs), idxs, dists, missed
end

