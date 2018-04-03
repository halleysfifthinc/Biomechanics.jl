export Trial,
       RawTrial,
       AnalyzedTrial,
       TrialDescriptor

export readtrial

"""
    TrialDescriptor

In recognition that data sources vary, and the mechanism of reading trials will differ 
between data sources, implement a TrialDescriptor for your data.
"""
abstract type TrialDescriptor end

"""
    Trial{TD}

A trial should describe the particulars of the data source. Trials are parameterized for
different datasources to allow for dispatching by the Trial parameter.
"""
struct Trial{TD}
    "The subject identifier"
    subject::Integer
    "The trial name"
    name::String
    "The absolute path to the file containing the trial data"
    path::String
    "The specific trial conditions; if unneeded, this can be empty"
    conds::Dict{Symbol,Symbol}

    """
        Trial{TD}(path[, subbase])

    Create a trial and infer the subject id and trialname from the path. Create an empty set
    of conditions.
    """
    function Trial{TD}(path::String,
                   subbase::String="Subject") where TD <: TrialDescriptor
        isabspath(path) || throw(ArgumentError("path must be absolute"))

        name = splitext(basename(path))[1]

        m = match(Regex("[\\\\,\\/]($(strip(subbase)) (?<subject>\\d+))"), path)
        m == nothing && throw(DomainError("no matching subject ID found for the given subbase and path"))

        return new(parse(m[:subject]), name, path, Dict{Symbol,Symbol}())
    end

    """
        Trial{TD}(subject, name, path[, conds])

    Create a completely specified trial.
    """
    function Trial{TD}(s,n,p,conds=Dict{Symbol,Symbol}()) where TD <: TrialDescriptor
        isabspath(p) || throw(ArgumentError("path must be absolute"))

        return new(s,n,p,conds)
    end
end

Base.show(io::IO, t::Trial) = print(io, t.subject, ", ", t.name, ", ", t.conds)

function Base.show(io::IO, ::MIME"text/plain", t::Trial{TD}) where TD
    println(io, "Trial{",TD,"}")
    println(io, "Subject => ", t.subject)
    println(io, "Name    => ", t.name)
    println(io, "Conditions:\n  ", t.conds)
end

"""
    readtrial(::Trial{TD})

Throws a MethodError. `readtrial` methods must be defined for a particular TrialDescriptor.
"""
function readtrial(t::Trial{TD}) where TD
    throw(MethodError(readtrial, (t,)))
end

"""
    RawTrial

A RawTrial wraps the base Trial, and includes any data and events which may be of interest.
"""
struct RawTrial
    t::Trial
    events::Dict{Symbol,Array}
    data::Matrix
end

Base.show(io::IO, r::RawTrial) = print(io, '(', r.t, "), ", keys(r.events), ", ",
                                       typeof(r.data), size(r.data))

function Base.show(io::IO, ::MIME"text/plain", r::RawTrial)
    show(io, MIME("text/plain"), r.t)
    println(io, "Events:")
    for key in keys(r.events)
        println(io, "  ", key)
    end
    println(io, "Data:")
    println(io, "  ", typeof(r.data), size(r.data))
end

"""
    AnalyzedTrial

Contains the results of any analysis/analyses performed on the trial.
"""
struct AnalyzedTrial
    t::Trial
    results::Dict{Symbol,Any}
end
