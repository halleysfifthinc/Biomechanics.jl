module Biomechanics

using DSP, Interpolations, Statistics

export timenormalize, normalizeeevents, normtime, avgcycle, demean, demean!, detrend,
       detrend!, centraldiff, increasing, decreasing, circmean, circstd, xcom,
       slidingwindow, calcresiduals, optfc

export ForwardBackwardPad

include("timenormalize.jl")
include("helperfuncs.jl")

end # module
