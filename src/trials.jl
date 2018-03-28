export Trial,
       RawTrial,
       AnalyzedTrial,
       TrialDescriptor

export readtrial

abstract type TrialDescriptor end

struct Trial{TD}
    subject::Int
    name::String
    path::String
    conds::Dict{Symbol,Symbol}

    function Trial{TD}(path::String;
                   subbase::String="Subject") where TD <: TrialDescriptor
        isabspath(path) || throw(ArgumentError("path must be absolute"))

        name = splitext(basename(path))[1]

        m = match(Regex("[\\\\,\\/]($subbase (?<subject>\\d+))"), path)

        return new(parse(m[:subject]), name, path, Dict{Symbol,Symbol}())
    end

    function Trial{TD}(s,n,p,conds;
                   subbase::String="Subject") where TD <: TrialDescriptor
        isabspath(p) || throw(ArgumentError("path must be absolute"))
        @assert n == splitext(basename(p))[1]

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

function readtrial(::Trial{TD}) where TD
    throw(MethodError("no matching method for `readtrial(::Trial{$TD})`. "*
                      "Extend `readtrial` to add support for you TrialDescriptor."))
end

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

struct AnalyzedTrial
    t::Trial
    results::Dict{Symbol,Any}
end
