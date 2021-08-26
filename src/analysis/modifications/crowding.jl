
"""
    add_crowding_constraint(weights::Dict{Int64, Float64})

Adds a molecular crowding constraint to the optimization problem: `∑ wᵢ × vᵢ ≤ 1` where `wᵢ`
is a weight and `vᵢ` is a flux index in the model's reactions specified in `weights` as `vᵢ
=> wᵢ` pairs.

See Beg, Qasim K., et al. "Intracellular crowding defines the mode and sequence of
substrate uptake by Escherichia coli and constrains its metabolic activity." Proceedings of
the National Academy of Sciences 104.31 (2007) for more details.
"""
add_crowding_constraint(weights::Dict{Int64,Float64}) =
    (model, opt_model) -> begin
        idxs = collect(keys(weights)) # order of keys and values is the same
        ws = values(weights)
        # since fluxes can be positive or negative, need absolute value: ∑ wᵢ × |vᵢ| ≤ 1
        # introduce slack variables to handle this
        @variable(opt_model, crowding_slack[1:length(weights)])
        @constraint(opt_model, crowding_slack .>= opt_model[:x][idxs])
        @constraint(opt_model, crowding_slack .>= -opt_model[:x][idxs])
        @constraint(opt_model, sum(w * crowding_slack[i] for (i, w) in enumerate(ws)) <= 1)
    end

"""
    add_crowding_constraint(weight::Float64; kwargs)

Variant of [`add_crowding_constraint`](@ref) that takes a single weight and assigns it to
each internal reaction flux, where internal reactions are identified with
[`find_internal_reactions`](@ref) and `kwargs` are passed to this function.
"""
add_crowding_constraint(weight::Float64; kwargs...) =
    (model, opt_model) -> begin
        idxs = find_internal_reactions(model; kwargs...)
        add_crowding_constraint(Dict(zip(idxs, fill(weight, length(idxs)))))(
            model,
            opt_model,
        )
    end

"""
    add_crowding_constraint(weights::Dict{String, Float64})

Variant of [`add_crowding_constraint`](@ref) that takes a dictinary of reactions `ids`
instead of reaction indices mapped to weights.
"""
add_crowding_constraint(weights::Dict{String,Float64}) =
    (model, opt_model) -> begin
        idxs = indexin(keys(weights), reactions(model))
        nothing in idxs && throw(ArgumentError("Reaction id not found in model."))
        add_crowding_constraint(Dict(zip(Int.(idxs), values(weights))))(model, opt_model)
    end