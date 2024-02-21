module Biomechanics

using Reexport, DSP, Interpolations, Statistics, StatsBase, LinearAlgebra, NearestNeighbors,
    Distances, LoopVectorization, Peaks

# utils.jl
export intervals

#reductions.jl
export circmean, circmeand, circstd, circstdd, mean_std_range, avgextrema, intervalextrema

# transformations.jl
export continuousphase, relativephase, crpensemble, demean, demean!, detrend, detrend!,
    centraldiff
export ForwardBackwardPad

# timenormalize.jl
export timenormalize, timestoindices, normtime, limitcycle, timeintervals_toindices,
    time_denormalize

# gaitevents.jl
export matchevents, toindices, totimes

# gaitstability.jl
export xcom, margin_of_stability, phase_coordination_index

# helperfuncs.jl
export calcresiduals, optfc

include("utils.jl")
include("reductions.jl")
include("transformations.jl")
include("timenormalize.jl")

include("stepmetrics.jl")
@reexport using .SpatiotemporalMetrics

include("gaitevents.jl")
include("gaitstability.jl")
include("helperfuncs.jl")

end # module
