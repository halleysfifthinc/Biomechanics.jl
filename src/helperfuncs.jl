using DSP

export calcresiduals,
       optfc,
       demean,
       demean!,
       detrend,
       detrend!,
       increasing,
       decreasing,
       slidingwindow

"""
    demean(x)

Subtract the mean of `x` from `x`
"""
function demean(x)
    return x .- mean(x)
end

"""
    demean!(x)

Mutate `x` by subtracting its mean
"""
function demean!(x)
    x .-= mean(x)
end

"""
    detrend(y)

Remove the linear trend from `y`
"""
function detrend(y)
    a, b = linreg(1:length(y),y)
    @. return y - (a + b*y)
end

"""
    detrend!(y)

Mutate `y` by removing its linear trend
"""
function detrend!(y)
    a, b = linreg(1:length(y),y)
    @. y -= (a + b*y)
end

"""
    decreasing(x)

Check if a vector is monotonically decreasing
"""
function decreasing(x)
    i = iterate(x)
    i === nothing && return false
    (last, state) = i
    i = iterate(x, state)

    while i !== nothing
        (next, state) = i
        next <= last || return false

        i = iterate(x, state)
        last = next
    end

    return true
end

"""
    increasing(x)

Check if a vector is monotonically increasing
"""
function increasing(x)
    i = iterate(x)
    i === nothing && return false
    (last, state) = i
    i = iterate(x, state)

    while i !== nothing
        (next, state) = i
        next >= last || return false

        i = iterate(x, state)
        last = next
    end

    return true
end

"""
    slidingwindow(f, A, len)

Apply function `f` to a sliding window of length `len` along `A`
"""
function slidingwindow(f, A::AbstractArray{T}, len) where T
    @assert length(A) >= len

    top = length(A) - len
    res = Array{T}(undef, top + 1)
    rg = 1:len

    for i in 0:top
        res[i+1] = f(@view(A[rg .+ i]))
    end

    return res
end

"""
    calcresiduals(data::Vector{Float64},
                    fs::Int64,
                    fc::Vector{Float64}=collect(.05:.05:25.0))::Vector{Float64}

Calculate the residuals of the data for a range of cutoff frequencies.

Following the method outlined by Winter (1990), the residual for a given cutoff freqency is defined as:
```math
R(f_c) = \\sqrt{\\frac{1}{N}\\sum_{i=1}^{N}(X_i - X̂_i)^2}
```
where ``f_c`` is the cutoff frequency, ``X_i`` is raw data at the *i* th sample, and ``X̂_i`` is filtered data at the *i* th sample [^1].

[1]D.A. Winter, Choice of Cutoff Frequency - Residual Analysis, in: Biomechanics and Motor Control of Human Movement, 4th ed., John Wiley & Sons, Inc., Hoboken, New Jersey, 2009: pp. 70–73.
"""
function calcresiduals(data::Vector{Float64},
                         fs::Real,
                         fc::Vector{Float64}=collect(.05:.05:25.0))::Vector{Float64}
    designmethod = Butterworth(2)
    len = length(fc)
    R = zeros(Float64,len)
    N = length(data)
    for j in 1:length(fc)
        responsetype = Lowpass(fc[j]*1.247; fs=fs)
        thefilter = digitalfilter(responsetype, designmethod)
        R[j] = sqrt((1/N)*sum((data.-filtfilt(thefilter,data)).^2))
    end
    return R
end

"""
    optfc(R::Vector{Float64},
                    fc::Vector{Float64}=collect(.05:.05:25.0),
                    lb::Float64=15.0)::Tuple{Float64,Float64,Float64}

Choose the optimal cutoff frequency based on the residuals.

The optimal cutoff frequency is the point at which the signal distortion is equal to the residual noise. This can be calculated by calculating the intercept at 0Hz for random noise, and then using the frequency at which the residuals fall below this intercept as the cutoff frequency [^1].

[1]D.A. Winter, Choice of Cutoff Frequency - Residual Analysis, in: Biomechanics and Motor Control of Human Movement, 4th ed., John Wiley & Sons, Inc., Hoboken, New Jersey, 2009: pp. 70–73.

"""
function optfc(R::Vector{Float64},
                         fc::Vector{Float64}=collect(.05:.05:25.0),
                         lb::Float64=15.0)::Tuple{Float64,Float64,Float64}
    q = findin(fc,lb)
    c = q[1]
    a,b = linreg(fc[c:end],R[c:end])
    i = findfirst(x -> x <= a, R)
    if i == 0
        return (0.0, a, b)
    else
        return (fc[i], a, b)
    end
end

