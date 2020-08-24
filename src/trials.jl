using Glob
using Glob: GlobMatch

export Trial, Segment, AnalyzedSegment, AbstractDataSource, DataFileType, TrialConditions

export readtrial, findtrials

struct DatasetSpec{DS<:AbstractDataSource}
    name::String
    datasource::DS
    pattern::String
end

function DatasetSpec(name, datasource::DS, pattern) where DS
    return DatasetSpec{DS}(name, datasource, pattern)
end

struct TrialConditions
    condnames::Vector{Symbol}
    required::Vector{Symbol}
    labels_rg::Regex
    subst::Vector{Pair{Regex,String}}
    types::Vector{Type}
end

"""
    TrialConditions(conditions, labels; kwargs...) -> TrialConditions

- `conditions` is a collection of condition names (eg `(:medication, :strength)`)
- `labels` is a `Dict` with keys for each condition name (eg `haskey(labels, :medication)`). Each key gets a collection of the labels for all levels and any transformation desired for that condition.

# Keyword arguments

- `required`: The conditions which every trial must have (in the case of some trials having optional/additional conditions)
- `types`: The (Julia) types for each condition (eg `[String, Int]`)
- `sep`: The character separating condition labels
"""
function TrialConditions(
    conditions,
    labels;
    required=conditions,
    types=fill(String, length(conditions)),
    sep="[_-]"
)
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

    return TrialConditions(collect(conditions), collect(required), Regex(labels_rg), subst, types)
end

"""
    Trial{S}

A `Trial` describes the referenced trial. Trials are parameterized for
different locations to allow for dispatching by the Trial parameter.
"""
struct Trial{S}
    "The subject identifier"
    subject::S

    "The trial name"
    name::String

    "The paths to any associated files containing trial data"
    paths::Dict{String,String}

    "The source type of the `paths`"
    sources::Dict{String,<:AbstractSource}

    "The specific trial conditions; if unneeded, this can be empty"
    conds::Dict{Symbol}
end

function Trial(subject::S, name, paths, sources, conds) where S
    return Trial{S}(subject, name, paths, sources, conds)
end

function findtrials(
    datasets::AbstractVector{DatasetSpec},
    conditions::TrialConditions;
    sid_type::Type=Int,
    subject_fmt=r"(?:Subject (?<subject>\d+))",
    ignorefiles::Union{Nothing, Vector{String}}=nothing,
    defaultconds::Union{Nothing, Dict{Symbol}}=nothing
)
    trials = Vector{Trial{sid_type}}()
    rg = subject_fmt*r".*"*conditions.labels_rg
    reqcondnames = conditions.required
    optcondnames = setdiff(conditions.condnames, reqcondnames)
    _defaultconds = Dict(cond => nothing for cond in conditions.condnames)
    if !isnothing(defaultconds)
        merge!(_defaultconds, defaultconds)
    end

    for set in datasets
        pattern = set.pattern
        files = glob(pattern)
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
                    all(enumerate(conditions.condnames)) do (i, cond)
                        trialcond = get(trial.conds, cond, _defaultconds[cond])
                        if isnothing(m[cond])
                            return isnothing(trialcond)
                        elseif conditions.types[i] === String
                            return m[cond] == trialcond
                        else
                            parse(conditions.types[i], m[cond]) == trialcond
                        end
                    end
                end

                if isempty(seenall)
                    conds = Dict(cond => String(m[cond]) for cond in reqcondnames)
                    foreach(enumerate(optcondnames)) do (i, cond)
                        if !isnothing(m[cond])
                            if conditions.types[i] === String
                                conds[cond] = String(m[cond])
                            else
                                conds[cond] = parse(conditions.types[i], m[cond])
                            end
                        end
                    end
                    push!(trials, Trial(sid, name, Dict(set.name => file),
                        Dict(set.name => set.datasource), conds))
                else
                    seen = only(seenall)
                    t = trials[seen]
                    if haskey(t.paths, set.name)
                        # println(repr(t.paths[set.first]))
                        # TODO: Implement `DuplicateTrialSourceError` and informative error msg
                        @show t.paths[set.name] file
                        continue
                    else
                        t.paths[set.name] = file
                        t.sources[set.name] = set.datasource
                    end
                end
            end
        end
    end

    return trials
end

function Base.show(io::IO, t::Trial)
    print(io, "Trial(", repr(t.subject), ", ", repr(t.name),
        ", $(length(t.sources)) sources, ", t.conds, ')')
end

function Base.show(io::IO, ::MIME"text/plain", t::Trial{S}) where D
    println(io, "Trial{", S, "}")
    println(io, "  Subject: ", t.subject)
    println(io, "  Name: ", t.name)
    print(io, "  Sources:")
    for p in t.sources
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
struct Segment{S,ID}
    trial::Trial{ID}
    conds::Dict{Symbol}
    source::S
    stime::Float64 # Start time
    etime::Float64 # End time
    metadata::Dict{String} # Is there a better key type or container? Symbol?

    function Segment{S}(trial::Trial{S}, conds::Dict{Symbol}, data) where S
        return new(trial, merge(conds, trial.conds), data)
    end
end

function Segment(trial::Trial{S}, conds::Dict{Symbol}, data) where S
    return Segment{S}(trial, conds, data)
end

function Base.show(io::IO, s::Segment{DS}) where DS
    print(io, "Segment{",DS,"}(",s.trial,",", s.conds, ")")
end

function Base.show(io::IO, mimet::MIME"text/plain", s::Segment{DS}) where DS
    print(io, "Segment{",DS,"}\n  ")
    show(io, s.trial)
    show(io, mimet, s.conds)
    show(io, mimet, s.data)
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

