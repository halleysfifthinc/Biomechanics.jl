using Glob
using Glob: GlobMatch

export Trial, Segment, AnalyzedSegment, AbstractDataSource, DataFileType, TrialConditions

export readtrial, findtrials

"""
    AbstractDataSource

In recognition that data sources vary, and the mechanism of reading trials will differ
between data sources, implement a subtype of AbstractDataSource for your data.
"""
abstract type AbstractDataSource end

struct DataFileType{R,P}
    rootdir::R
    gpath::P
end

struct TrialConditions
    condnames::Vector{Symbol}
    required::Vector{Symbol}
    labels_rg::Regex
    subst::Vector{Pair{Regex,String}}
end

function TrialConditions(conditions, labels; required=conditions, sep="[_-]")
    labels_rg = ""
    subst = Vector{Pair{Regex,String}}(undef, 0)

    for cond in conditions
        labels_rg *= "(?<$cond>"
        if labels[cond] isa Regex
            labels_rg *= labels[cond].pattern
        else
            labels_rg *= join((x isa Pair ? x.second : x for x in labels[cond]), '|')
        end
        optchar = cond in required ? "" : "?"
        labels_rg *= string(')', optchar, sep, '?')
        foreach(labels[cond]) do condlabel
            if condlabel isa Pair
                altlabels = condlabel.first isa Union{Symbol,String} ? [condlabel.first] : condlabel.first
                filter!(label -> label != condlabel.second, altlabels)
                push!(subst, Regex("(?:"*join(altlabels, '|')*")") => condlabel.second)
            end
        end
    end

    return TrialConditions(collect(conditions), collect(required), Regex(labels_rg), subst)
end

"""
    Trial{DS}

A `Trial` describes the referenced trial. Trials are parameterized for
different datasources to allow for dispatching by the Trial parameter.
"""
struct Trial{DS<:AbstractDataSource,S}
    "The subject identifier"
    subject::S

    "The trial name"
    name::String

    "The absolute paths to any associated files containing trial data"
    paths::Dict{String,String}

    "The specific trial conditions; if unneeded, this can be empty"
    conds::Dict{Symbol}
end

function findtrials(
    DS,
    datasources,
    conditions;
    sid_type::Type=Int,
    subject_fmt=r"(?:Subject (?<subject>\d+))",
    ignorefiles::Union{Nothing, Vector{String}}=nothing,
    defaultconds::Union{Nothing, Dict{Symbol}}=nothing
)
    trials = Vector{Trial{DS,sid_type}}()
    rg = subject_fmt*r".*"*conditions.labels_rg
    reqcondnames = conditions.required
    optcondnames = setdiff(conditions.condnames, reqcondnames)
    _defaultconds = Dict(cond => nothing for cond in conditions.condnames)
    if !isnothing(defaultconds)
        merge!(_defaultconds, defaultconds)
    end

    for ds in datasources
        rootdir = ds.second.rootdir
        gpath = ds.second.gpath
        files = glob(gpath, rootdir)
        if !isnothing(ignorefiles)
            setdiff!(files, ignorefiles)
        end

        for file in files
            _file = foldl((str, pat) -> replace(str, pat), conditions.subst; init=file)
            m = match(rg, _file)
            isnothing(m) && continue

            if isnothing(m[:subject]) || any(isnothing.(m[cond] for cond in reqcondnames))
                continue
            else
                name = splitext(basename(file))[1]
                sid = !(sid_type <: String) ? parse(sid_type, m[:subject]) : String(m[:subject])
                seenall = findall(trials) do trial
                    trial.subject == sid &&
                    all(m[cond] == get(trial.conds, cond, _defaultconds[cond])
                        for cond in conditions.condnames)
                end

                if isempty(seenall)
                    conds = Dict(cond => String(m[cond]) for cond in reqcondnames)
                    foreach(optcondnames) do cond
                        if !isnothing(m[cond])
                            conds[cond] = String(m[cond])
                        end
                    end
                    push!(trials, Trial{DS,sid_type}(sid, name, Dict(ds.first => file),
                        conds))
                else
                    seen = only(seenall)
                    t = trials[seen]
                    if haskey(t.paths, ds.first)
                        # println(repr(t.paths[ds.first]))
                        # TODO: Implement `DuplicateTrialSourceError` and informative error msg
                        @show t.paths[ds.first] file
                        continue
                    else
                        t.paths[ds.first] = file
                    end
                end
            end
        end
    end

    return trials
end

function Base.show(io::IO, t::Trial)
    print(io, "Trial(", repr(t.subject), ", ", repr(t.name), ", $(length(t.paths)) paths, ", t.conds, ')')
end

function Base.show(io::IO, ::MIME"text/plain", t::Trial{DS,S}) where {DS,S}
    println(io, "Trial{", DS, ',', S, "}")
    println(io, "  Subject: ", t.subject)
    println(io, "  Name: ", t.name)
    print(io, "  Paths:")
    for p in t.paths
        print(io, "\n    ")
        show(io, p)
    end
    print(io, "\n  Conditions:")
    for c in t.conds
        print(io, "\n    ")
        show(io, c)
    end
    println(io)
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

Base.show(io::IO, s::Segment{DS}) where DS = print(io, "Segment{",DS,"}(",s.trial,",", s.conds, ")")

function Base.show(io::IO, ::MIME"text/plain", s::Segment{DS}) where DS
    print(io, "Segment{",DS,"}\n  ")
    show(io, s.trial)
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
    show(io, as.s)
    show(io, MIME("text/plain"), as.results)
end

