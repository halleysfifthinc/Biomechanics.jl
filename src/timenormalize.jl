using Interpolations

export timenormalizestrides,
       limitcycle


function timenormalizestrides(trial::Trial, st::Float64, en::T=nothing; kwargs...) where {T<:Void}
    res, rhs = innertimenormalizestrides(trial, st; kwargs...)

    return res   
end

function timenormalizestrides(trial::Trial, st::Float64, en::T; kwargs...) where T
    res, rhs = innertimenormalizestrides(trial, st; kwargs...)
    
    tststr = normtime(rhs[1],rhs[2])
    pst = findfirst(x -> x >= st*fs, tststr)
    pen = 0
    for s in 1:(length(rhs)-1)
        tststr = normtime(rhs[s],rhs[s+1])
        pen = findfirst(x -> x >= en*fs, tststr)+(s-1)*fs
        !iszero(mod(pen,100)) && break
    end
    
    return (res, pst, pen)
end

function innertimenormalizestrides(trial::Trial, st::Float64; kwargs...)

    rt = readtrial(trial, st; kwargs...)
    prerhs = rt.events[:RHS][1]
    lastrhs = rt.events[:RHS][end]
    cols = size(rt.data,2)

    res = zeros((lastrhs-prerhs)*100, cols)

    fill_normdims!(res, rt.data, rt.events[:RHS])

    return (res, rt.events[:RHS])
end

function isnanbyrows(x::Matrix)::UnitRange{Int64}
    cols = 1:size(x,2)
    sigstart = maximum(( findnext(!isnan, view(x, :,i), 1) for i in cols ))
    sigend = minimum(( findprev(!isnan, view(x, :,i), size(x,1)) for i in cols ))

    return sigstart:sigend
end

function fill_normdims!(res::AbstractArray, 
                        data::AbstractArray,
                        events::Vector{Int})
    (size(res,2) != size(data,2)) && throw(
        ArgumentError("`res` and `data` must have the same number of columns"))
    for i in 1:size(data,2)
        # Create interpolation object
        itp = interpolate(data[:,cols[i]], BSpline(Cubic(Line())), OnGrid())

        fill_normstrides!(view(res, :, i), itp, events)
    end

    nothing
end

function fill_normstrides!(nstr::AbstractArray, str::AbstractArray, events::Vector)
    @assert mod(size(nstr,1),100) === 0
    @assert size(nstr,1) === (length(events)-1)*100

    # Create time normalized strides
    @inbounds for s in 1:length(events)-1
        nstr[(1:100)+(s-1)*100] .= str[normtime(events[s],events[s+1])]
    end

    nothing
end

normtime(t1,t2) = linspace(t1,t2,101)[1:100]

"""
    limitcycle(rootdir::String, sub::Int)

Calculate the ensemble average of the exported data of the subject's steady-state trial

Returns the average and standard deviation time normalized as 100% the length of a stride.
"""
function limitcycle(trial::Trial, numstrides=100, st=30.0; kwargs...)
    nstr = timenormalizestrides(trial, st; numstrides=numstrides, kwargs...)

    NW = Array{Float64,2}(100, size(nstr,2))
    vNW = Array{Float64,2}(100, size(nstr,2))

    for dim in 1:size(nstr,2)
        NW[:,dim] .= [ mean(nstr[i:100:end, dim]) for i in 1:100 ]
        vNW[:,dim] .= [ std(nstr[i:100:end, dim]) for i in 1:100 ]
    end

    return (NW, vNW)
end