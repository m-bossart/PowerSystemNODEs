# TODO: Use Dict temporarily during dev while the fields are defined
struct NODETrainInputs
    inputs_name::String
    data::Dict{Symbol, Vector{Float64}}
end

function NODETrainInputs(name::String)
    return NODETrainInputs(name, Dict{Symbol, Vector{Float64}}())
end

function serialize(inputs::NODETrainInputs, file_path::AbstractString)
    open(file_path, "w") do io
        JSON3.write(io, inputs)
    end
    return
end
