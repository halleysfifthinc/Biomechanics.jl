module Biomechanics

using Reexport, DSP, Interpolations, Statistics, LinearAlgebra

# utils.jl
export intervals

#reductions.jl
export circmean, circstd, mean_std_range, avgextrema, intervalextrema

# transformations.jl
export continuousphase, demean, demean!, detrend, detrend!, centraldiff

# timenormalize.jl
export timenormalize, timestoindices, normtime, limitcycle


# gaitstability.jl
export xcom

# helperfuncs.jl
export calcresiduals, optfc

export ForwardBackwardPad

include("utils.jl")
include("reductions.jl")
include("transformations.jl")
include("timenormalize.jl")

include("stepmetrics.jl")
@reexport using .SpatiotemporalMetrics

include("gaitstability.jl")
include("helperfuncs.jl")

end # module
