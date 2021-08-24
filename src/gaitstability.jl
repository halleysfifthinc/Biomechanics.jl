"""
    xcom(com, vel, marker)

Compute the extrapolated center of mass.

`com`, and `marker` are all assumed to be 3-dimensional variables where the order of
dimensions is `[ML, AP, VT]`. Units are assumed to be standard SI units.

See McAndrew et al. (2012)[https://doi.org/10.1016/j.jbiomech.2011.12.027]
"""
function xcom(com, bos)
    size(com, 1) === size(bos, 1) || throw(ArgumentError("Both signals must have equal length"))

    # tsum = sum((com - bos).^2; dims=2)
    # l = mean(sqrt, tsum)
    d = com - bos
    l = mean(filter(!isnan, mapslices(norm, d; dims=2)))
    ω₀ = sqrt(9.81/l)

    vel = centraldiff(com; padding=ForwardBackwardPad())

    xCOM = com + vel ./ ω₀

    return xCOM
end

function margin_of_stability(com, lbos, rbos, levents, revents; axis=1)
    lmos = lbos - xcom(com, lbos)
    rmos = rbos - xcom(com, rbos)
    sp = sortperm([levents; revents])

    lmos_ev, rmos_ev = lmos[levents, axis], rmos[revents, axis]
    return lmos_ev, rmos_ev, [lmos_ev; rmos_ev][sp]
end

function phase_coordination_index(;lfs, rfs, lfo, rfo)
    lswing = mean(swing(lfs, lfo))
    rswing = mean(swing(rfs, rfo))

    if lswing < rswing
        tₗ, tₛ = beginwithevent(rfs, lfs)
    else
        tₗ, tₛ = beginwithevent(lfs, rfs)
    end
    resize!(tₗ, min(length(tₗ), length(tₛ)))
    resize!(tₛ, min(length(tₗ), length(tₛ)))

    φ = 360 .* (tₛ[1:end-1] .- tₗ[1:end-1])/stridetimes(tₗ)
    φ_ABS = mean(abs, φ .- 180)
    Pφ_ABS = (φ_ABS/180)*100
    φ_CV = variation(φ)*100

    return φ_CV + Pφ_ABS
end

