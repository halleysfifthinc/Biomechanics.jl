module Biomechanics

using DSP, Interpolations, Statistics, LinearAlgebra

export timenormalize, timestoindices, normtime, limitcycle, demean, demean!, detrend,
       detrend!, centraldiff, circmean, circstd, xcom, calcresiduals, optfc

export ForwardBackwardPad

include("timenormalize.jl")
include("helperfuncs.jl")
include("stepmetrics.jl")

end # module
