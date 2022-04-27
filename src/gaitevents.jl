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

function roerdink2008(; lheel, rheel, fs, fc_minprom=30, fo_minprom=1200)
    lfcpred, _ = peakproms!(argminima(lheel, 10), lheel; minprom=fc_minprom)
    rfcpred, _ = peakproms!(argminima(rheel, 10), rheel; minprom=fc_minprom)

    lheel_vel = centraldiff(lheel; dt=inv(fs), padding=ForwardBackwardPad())
    rheel_vel = centraldiff(rheel; dt=inv(fs), padding=ForwardBackwardPad())

    lfopred, _ = peakproms!(argmaxima(lheel_vel, 10), lheel_vel; minprom=fo_minprom)
    rfopred, _ = peakproms!(argmaxima(rheel_vel, 10), rheel_vel; minprom=fo_minprom)

    return (;lfs=totimes(lfcpred, fs), rfs=totimes(rfcpred, fs),
             lfo=totimes(lfopred, fs), rfo=totimes(rfopred, fs))
end

