"""
    xcom(pos, vel, marker)

Compute the extrapolated center of mass.

`pos`, `vel`, and `marker` are all assumed to be 3-dimensional variables where the order of
dimensions is `[ML, AP, VT]`. Units are assumed to be standard SI units.

See McAndrew et al. (2012)[https://doi.org/10.1016/j.jbiomech.2011.12.027]
"""
function xcom(pos, vel, marker)
    size(pos, 1) === size(vel, 1) || throw(ArgumentError("Both signals must have equal length"))

    tsum = sum((pos - marker).^2; dims=2)
    l = mean(sqrt.(tsum))
    ω₀ = sqrt(9.81/l)

    ML = pos[:,1] + vel[:,1]./ω₀
    AP = pos[:,2] + vel[:,2]./ω₀

    return (AP, ML)
end

