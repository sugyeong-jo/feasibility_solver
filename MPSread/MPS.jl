__precompile__()

module MPS

# order, indices, weights
const SOS = Tuple{Int, Vector{Int}, Vector{Float64}}

const START_REG = r"^([\.0-9eE])"
const NAME_REG = r"([^a-zA-Z0-9\!\"\#\$\%\&\(\)\/\,\.\;\?\@\_\`\'\{\}\|\~])"

function verifyname(name::String)
    if length(name) > 16
        return false
    end
    m = match(START_REG, name)
    if !isa(m, Void)
        return false
    end
    m = match(NAME_REG, name)
    if !isa(m, Void)
        return false
    end
    return true
end

function correctname(name::String)
    if length(name) > 16
        warn("Name $(name) too long. Truncating.")
        return correctname(String(name[1:16]))
    end
    m = match(START_REG, name)
    if !isa(m, Void)
        warn("Name $(name) cannot start with a period, a number, e, or E. Removing from name.")
        return correctname(replace(name, START_REG, ""))
    end
    m = match(NAME_REG, name)
    if !isa(m, Void)
        warn("Name $(name) contains an illegal character. Removing from name.")
        return correctname(replace(name, NAME_REG, ""))
    end
    return name
end

include("writer.jl")
include("reader.jl")

end