using Interpolations

export timenormalize,
       avgcycle


function timenormalize(data::AbstractArray, events::AbstractVector,len=100)
    cols = size(data,2)
    res = zeros((length(events)-1)*len, cols)

    fill_normdims!(res, data, events; len=len)

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
    avgcycle(data, len=100)

Calculate the ensemble average and stdev of `data` assuming a normalized length of `len`
"""
function avgcycle(data::AbstractArray, len=100)
    mod(size(data,1),len) == 0 || throw(ArgumentError("length of data must be even"*
                                                      " multiple of len"))
    normal = Array{Float64,2}(len, size(nstr,2))
    vnormal = Array{Float64,2}(len, size(nstr,2))

    for dim in 1:size(nstr,2)
        normal[:,dim] .= [ mean(nstr[i:len:end, dim]) for i in 1:len ]
        vnormal[:,dim] .= [ std(nstr[i:len:end, dim]) for i in 1:len ]
    end

    return (normal, vnormal)
end
