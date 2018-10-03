export Trial,
       Segment,
       AnalyzedSegment,
       AbstractDataSource

export readtrial

"""
    AbstractDataSource

In recognition that data sources vary, and the mechanism of reading trials will differ 
between data sources, implement a subtype of AbstractDataSource for your data.
"""
abstract type AbstractDataSource end

"""
    Trial{DS}

A `Trial` describes the referenced trial. Trials are parameterized for
different datasources to allow for dispatching by the Trial parameter.
"""
struct Trial{DS<:AbstractDataSource}
    "The subject identifier"
    subject::Integer

    "The trial name"
    name::String

    "The absolute path to the file containing the trial data"
    path::String

    "The specific trial conditions; if unneeded, this can be empty"
    conds::Dict{Symbol}

    """
        Trial{DS}(path[, subbase])

    Create a trial and infer the subject id and trialname from the path. Create an empty set
    of conditions.
    """
    function Trial{DS}(path::String,
                   subbase::String="Subject") where DS <: AbstractDataSource
        isabspath(path) || throw(ArgumentError("path must be absolute"))

        name = splitext(basename(path))[1]

        m = match(Regex("[\\\\,\\/]($(strip(subbase)) (?<subject>\\d+))"), path)
        m == nothing && throw(DomainError("no matching subject ID found for the given subbase and path"))

        return new(tryparse(Int, m[:subject]), name, path, Dict{Symbol,Any}())
    end

    """
        Trial{DS}(subject, name, path[, conds])

    Create a completely specified trial.
    """
    function Trial{DS}(s,n,p,conds=Dict{Symbol,Any}()) where DS <: AbstractDataSource
        isabspath(p) || throw(ArgumentError("path must be absolute"))

        return new(s,n,p,conds)
    end
end

Base.show(io::IO, t::Trial) = print(io, t.subject, ", ", t.name, ", ", t.conds)

function Base.show(io::IO, ::MIME"text/plain", t::Trial{DS}) where DS
    println(io, "Trial{",DS,"}")
    println(io, "Subject => ", t.subject)
    println(io, "Name    => ", t.name)
    println(io, "Conditions:\n  ", t.conds)
end

function readtrial(t::Trial{DS}) where DS
    throw(MethodError(readtrial, (t,)))
end

"""
    Segment

A `Segment` is a container which includes all of, or a part of the data from a particular `Trial`. 
"""
struct Segment{DS<:AbstractDataSource}
    trial::Trial{DS}
    conds::Dict{Symbol}
    data::DS
end

Base.show(io::IO, s::Segment{DS}) where DS = print(io, "Segment{",DS,"}(",s.trial,",", s.conds, ",", s.data, ")")

function Base.show(io::IO, ::MIME"text/plain", s::Segment{DS}) where DS
    println(io, "Segment{",DS,"}")
    show(io, MIME("text/plain"), s.trial)
    show(io, MIME("text/plain"), s.conds)
    show(io, MIME("text/plain"), s.data)
end

"""
    AnalyzedSegment

Contains the results of any analysis/analyses performed on the trial segment.
"""
struct AnalyzedSegment{DS<:AbstractDataSource}
    s::Segment{DS}
    results::Dict{Symbol}
end

Base.show(io::IO, as::AnalyzedSegment) = print(io, "AnalyzedSegment(",as.s,",", as.results,")")

function Base.show(io::IO, ::MIME"text/plain", as::AnalyzedSegment{DS}) where DS
    println(io, "AnalyzedSegment{",DS,"}")
    show(io, MIME("text/plain"), as.s)
    show(io, MIME("text/plain"), as.results)
end

