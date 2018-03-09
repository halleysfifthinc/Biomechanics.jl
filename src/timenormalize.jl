using Interpolations

export timenormalizestrides,
       limitcycle


function timenormalizestrides(trial::Trial, st::Float64, en::T=nothing; kwargs...) where {T<:Void}
    rt = readtrial(trial, st; kwargs...)

    res = timenormalize(rt.data, rhs)

    return res   
end

function timenormalizestrides(trial::Trial, st::Float64, en::T; kwargs...) where T
    rt = readtrial(trial, st; kwargs...)
    rhs = rt.events[:RHS]

    res = timenormalize(rt.data, rhs)
    
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

function timenormalize(data::AbstractArray, events::AbstractVector,len=100)
    cols = size(data,2)
    res = zeros((length(events)-1)*len, cols)

    fill_normdims!(res, data, events; len=len)

    return res
end

function isnanbyrows(x::Matrix)::UnitRange{Int64}
    cols = 1:size(x,2)
    sigstart = maximum(( findnext(!isnan, view(x, :,i), 1) for i in cols ))
    sigend = minimum(( findprev(!isnan, view(x, :,i), size(x,1)) for i in cols ))

    return sigstart:sigend
end

function fill_normdims!(res::AbstractArray, 
                        data::AbstractArray,
                        events::Vector{Int},
                        len=100)
    (size(res,2) != size(data,2)) && throw(
        ArgumentError("`res` and `data` must have the same number of columns"))
    for i in 1:size(data,2)
        # Create interpolation object
        itp = interpolate(@view data[:,i], BSpline(Cubic(Line())), OnGrid())

        fill_normstrides!(@view res[:,i], itp, events; len=len)
    end

    nothing
end

function fill_normstrides!(nstr::AbstractArray, str::AbstractArray, events::Vector,len=100)
    @assert mod(size(nstr,1),len) === 0
    @assert size(nstr,1) === (length(events)-1)*len

    # Create time normalized strides
    @inbounds for s in 1:length(events)-1
        nstr[(1:len)+(s-1)*len] .= str[normtime(events[s],events[s+1],len=len)]
    end

    nothing
end

function normtime(t1,t2,len=100)
    st = (t2-t1)/len
    t1:st:(t2-st)
end

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