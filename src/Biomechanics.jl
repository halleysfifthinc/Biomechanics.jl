module Biomechanics

using DSP, Interpolations, Statistics, LinearAlgebra

export timenormalize, timestoindices, normtime, avgcycle, demean, demean!, detrend,
       detrend!, centraldiff, increasing, decreasing, circmean, circstd, xcom,
       slidingwindow, calcresiduals, optfc

export ForwardBackwardPad

include("timenormalize.jl")
include("helperfuncs.jl")
include("stepmetrics.jl")

end # module
