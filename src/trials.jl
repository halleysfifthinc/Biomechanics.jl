export Trial, RawTrial, AnalyzedTrial

export readtrial,
       getsessionorder

struct Trial
    subject::Int
    name::String
    path::String
    conds::Dict{Symbol,Symbol}

    function Trial(path::String;
                   subbase::String="Subject")
        isabspath(path) || throw(ArgumentError("path must be absolute"))
        ispath(path) || throw(ArgumentError("path must be existing file"))

        name = splitext(basename(path))[1]

        m = match(Regex("[\\\\,\\/]($subbase (?<subject>\\d\\d))"), path)

        return new(parse(m[:subject]), name, path, Dict{Symbol,Symbol}())
    end

    function Trial(s,n,p,conds;
                   subbase::String="Subject")
        isabspath(p) || throw(ArgumentError("path must be absolute"))
        ispath(p) || throw(ArgumentError("path must be existing file"))
        @assert n == splitext(basename(p))[1]

        return new(s,n,p,conds)
    end
end

Base.show(io::IO, t::Trial) = print(io, t.subject, ", ", t.name, ", ", t.conds)

function Base.show(io::IO, ::MIME"text/plain", t::Trial)
    println(io, "Subject => ", t.subject)
    println(io, "Name => ", t.subject)
    println(io, "Conditions:\n  ", t.conds)
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

const fs = 100

function readtrial(trial::Trial, st::Float64; kwargs...)
    kwargs == Dict(kwargs)

    if haskey(kwargs, :cols) && !isempty(kwargs[:cols])
        cols = kwargs[:cols]
    else
        cols = Vector{Int}()
        get(kwargs, :pos, true) && append!(cols, collect(2:4))
        get(kwargs, :linvel, true) && append!(cols, collect(5:7))
        get(kwargs, :ori, true) && append!(cols, collect(8:10))
        get(kwargs, :angvel, true) && append!(cols, collect(11:13))
    end

    data = readdlm(trial.path, '\t', Float64; skipstart=5)
    events = Dict{Symbol,Array}()
    rhs = round.(Int,filter(!isnan,trial[:,end])*fs)
    lhs = round.(Int,filter(!isnan,trial[:,end-1])*fs)

    prerhs = findfirst(x -> x >= st*fs,rhs)-1
    lastrhs = haskey(kwargs, :numstrides) ? (length(rhs)-1) : (prerhs+kwargs[:numstrides]::Int+1)

    try
        @assert lastrhs < length(rhs)
    catch
        print("\n")
        throw(ArgumentError("Insufficient number of strides in $sub, $tname"))
    end
    
    data = data[rhs[prerhs]:rhs[lastrhs],cols]
    events[:RHS] = rhs[prerhs:lastrhs]
    events[:LHS] = filter!(x -> (rhs[prerhs] <= x <= rhs[lastrhs]), lhs)

    return RawTrial(trial, events, data)
end

function getsessionorder(session::String)
    isabspath(session) || throw(ArgumentError("path must be absolute"))
    isdir(session) || throw(ArgumentError("path must be to a directory"))

    badenfs = Vector{Int}()
    trialctime = Vector{DateTime}()
    df = DateFormat("y,m,d,H,M,S")

    trials = readdir(session)
    filter!(file -> endswith(file,"Trial.enf"), trials)

    for trialnum in eachindex(trials)
        file = readstring(joinpath(session, trials[trialnum]))
        m = match(r"CREATIONDATEANDTIME\=(?<timestring>.*)\r", file)

        if m === nothing
            push!(badenfs, trialnum)
            continue
        end

        dt = DateTime(m[:timestring], df)
        push!(trialctime, dt)
    end

    deleteat!(trials, badenfs)

    order = sortperm(trialctime)

    return map(x -> basename(x[1:end-10]), trials[order])
end
