"""
    AbstractSource

In recognition that data sources vary, and the mechanism of reading trials will differ
between data sources, implement a subtype of AbstractDataSource for your data.
"""
abstract type AbstractSource end

function readsource end

struct C3DSource <: AbstractSource; end

function readsource(::C3DSource, filename::AbstractString; kwargs...)
    return readc3d(filename; kwargs...)
end

struct MATSource <: AbstractSource; end

function readsource(::MATSource,
    trial::Trial{<:AbstractDataSource},
    st::Float64 = 0.0,
    en::Float64 = Inf;
    events::Vector = [],
    ts::Vector = [],
    fs = 100
)
    isfile(trial.paths["export"]) || throw(ArgumentError("trial $(trial.paths["export"]) exist"))
    st >= 0.0 || throw(ArgumentError("start time must be positive"))
    st <= en || throw(ArgumentError("end time must be greater than start time"))

    evnames = Symbol.(events)
    file = matopen(trial.paths["export"])

    revents = Dict{Symbol,Vector{Float64}}()
    for e in events
        if exists(file, string(e))
            syme = Symbol(e)
            tmp = read(file, string(e))[1]
            if tmp isa AbstractArray
                revents[syme] = vec(tmp)
            else
                revents[syme] = [tmp]
            end
            strt = findfirst(x -> x >= st, revents[syme])
            if strt == 0
                # @warn "no $e events during given start and end times"
                delete!(revents, syme)
                break
            end
            endi = findlast(x -> x <= en, revents[syme])
            revents[syme] = revents[syme][strt:endi] .- st .+ (1/fs) # Shift events to be index accurate for data subsection
        else
            # @warn "Requested event $e does not exist in source data"
        end
    end

    data = Dict{Symbol,Matrix{Float64}}()
    for t in ts
        if exists(file, string(t))
            symt = Symbol(t)
            data[symt] = read(file, string(t))[1]
            len = size(data[symt], 1)
            strti = round(Int, st*fs)
            strti += (strti === 0) ? 1 : 0
            endi = en == Inf ? len : min(len, round(Int, en*fs))
            data[symt] = data[symt][strti:endi, :]
        else
            @warn "Requested time series $t does not exist in source data"
        end
    end

    return Segment(trial, Dict{Symbol,Any}(), DSData(revents, data))
end
