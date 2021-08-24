function rmsd_simd(a, b)
    N = length(a)
    N === length(b) || throw(DimensionMismatch())
    err = zero(promote_type(eltype(a), eltype(b)))
    @turbo for i ∈ eachindex(a,b)
        err += abs2(a[i] - b[i])
    end
    return sqrt(err/size(a,1))
end

logrange(x1, x2, n) = (exp10(y) for y in range(log10(x1), log10(x2), length=n))

"""
    calcresiduals(data, fc; fs)

Calculate the residuals of the data for a range of cutoff frequencies.

Following the method outlined by Winter (1990), the residual for a given cutoff freqency is defined as:

```math
R(f_c) = \\sqrt{\\frac{1}{N}\\sum_{i=1}^{N}(X_i - X̂_i)^2}
```
where ``f_c`` is the cutoff frequency, ``X_i`` is raw data at the *i* th sample, and ``X̂_i``
is filtered data at the *i* th sample [^1].

[1] D.A. Winter, Choice of Cutoff Frequency - Residual Analysis, in: Biomechanics and Motor
Control of Human Movement, 4th ed., John Wiley & Sons, Inc., Hoboken, New Jersey, 2009: pp.
70–73.

"""
function calcresiduals(data, fc=[range(1, fs/2/1.247, length=200);]; fs, db=false)
    sort!(fc)
    last(fc) ≤ fs/2/1.247 || throw(ArgumentError("cutoff frequency cannot be greater than Nyquist, corrected for dual-pass, at $(fs/2/1.247)"))
    numFc = length(fc)
    R = zeros(Float64,numFc)
    designmethod = Butterworth(2)
    filtdata = similar(data)

    for j in 1:numFc
        responsetype = Lowpass(fc[j]*1.247; fs=fs)
        thefilter = digitalfilter(responsetype, designmethod)
        filtdata[:] = filtfilt(thefilter, data)
        R[j] = rmsd_simd(data, filtdata)
    end

    if db
        baseP = rmsd_simd(data, zero(data))
        R[:] = amp2db.(R./baseP)
    end

    return R
end

"""
    optfc(data, fc; fs, [tol=0.9, rmsnoise])

Find the optimal cutoff frequency.

The optimal cutoff frequency is the point at which the signal distortion is equal to the
baseline RMS noise. The RMS noise is calculated as the intercept of the residual noise line,
which is the longest linear section of residuals, ending at the residual for the highest
frequency evaluated, with a correlation coefficient ≥`tol` [^1].

[1] D.A. Winter, Choice of Cutoff Frequency - Residual Analysis, in: Biomechanics and Motor
Control of Human Movement, 4th ed., John Wiley & Sons, Inc., Hoboken, New Jersey, 2009: pp.
70–73.
"""
function optfc(data, fc=[range(1, fs/2/1.247, length=200);]; fs, tol=0.9, rmsnoise=nothing)
    Rs = calcresiduals(data, fc; fs)
    if isnothing(rmsnoise)
        linregs = [ @views(((; a,b = _linreg(fc[n:end], Rs[n:end]))...,
                            r2 = cor((fc[n:end], Rs[n:end]))))
            for n in 1:length(fc)-1 ]
        d = findlast(x -> x.r2 < tol && x.b < 0, linregs)
        (isnothing(d) || d === lastindex(linregs)) &&
            throw(error("a linear section of residuals was not found with an r² ≥ $tol; try changing the tolerance"))
        rmsnoise = linregs[d+1].a
    end

    i = findlast(>(rmsnoise), Rs)
    @assert !isnothing(i) # At least one residual must be less than the intercept

    return fc[i], rmsnoise
end

