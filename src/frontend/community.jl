
# Copyright (c) 2021-2024, University of Luxembourg
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# TODO this might eventually deserve the QP version

# TODO steadycom

"""
$(TYPEDSIGNATURES)

Construct an instance of a linear community Flux Balance Analysis (cFBA) model.
The relative abundances of the organisms are known in advance; this function
predicts the maximum achievable community growth rate.

`model_abundances` is a dictionary-like object that maps model identifiers to
tuples of models (usually subtypes of `AbstractFBCModel`) and their abundances
(such as: `"bug1" => (bug1, 0.5), ...`). `community_exchange_bounds` is a
dictionary-like object that can be additionally used to restrict selected
community exchange reactions (keys should be reaction IDs, the values are
converted to `ConstraintTrees`-like bounds). Bounds otherwise default to
parameter `default_community_exchange_bound`, which itself defaults to
`nothing` (i.e., unbounded).

If required, constraint trees may be supplied instead of `AbstracFBCModel`s in
`model_abundances`. These must provide an interface compatible with
`interface_exchanges` and `interface_biomass`.

`interface` is forwarded to [`flux_balance_constraints`](@ref).
`interface_exchanges` and `interface_biomass` are used to pick up the correct
interface part to contribute to the community exchanges and community biomass.
"""
function community_flux_balance_constraints(
    model_abundances,
    community_exchange_bounds;
    interface = :identifier_prefixes,
    interface_exchanges = x -> x.interface.exchanges,
    interface_biomass = x -> x.interface.biomass,
    default_community_exchange_bound = nothing,
)
    @assert length(model_abundances) >= 1 "at least one community member is required"
    @assert isapprox(sum(a for (_, (_, a)) in model_abundances), 1) "community member abundances must sum to 1"

    bounds_lookup = Dict(community_exchange_bounds)

    constraints = interface_constraints(
        (
            Symbol(k) => if m isa A.AbstractFBCModel
                c = flux_balance_constraints(m; interface)
                (c, interface_exchanges(c), a)
            elseif m isa C.ConstraintTree
                (m, interface_exchanges(m), a)
            else
                throw(DomainError(m, "unsupported community member type"))
            end for (k, (m, a)) in model_abundances
        );
        out_interface = :community_exchanges,
        out_balance = :community_balance,
        bound = x ->
            get(bounds_lookup, String(last(x)), default_community_exchange_bound),
    )

    growth_sums = [
        Symbol(k) => C.Constraint(sum_value(interface_biomass(constraints[Symbol(k)])))
        for (k, _) in model_abundances
    ]

    constraints *
    :equal_growth^all_equal_constraints(
        last(growth_sums[1]).value,
        C.ConstraintTree(growth_sums[2:end]),
    ) *
    :community_biomass^last(growth_sums[1])
end

export community_flux_balance_constraints

"""
$(TYPEDSIGNATURES)

Run the Community Flux Balance Analysis on several models. All arguments are
forwarded to [`community_flux_balance_constraints`](@ref) which constructs the
model; this function returns the solved model.
"""
community_flux_balance_analysis(args...; kwargs...) = frontend_optimized_values(
    community_flux_balance_constraints,
    args...;
    objective = x -> x.community_biomass.value,
    kwargs...,
)

export community_flux_balance_analysis

"""
$(TYPEDSIGNATURES)

TODO TODO
"""
function community_composition_balance_constraints(
    models,
    growth,
    community_exchange_bounds;
    interface = :identifier_prefixes,
    interface_exchanges = x -> x.interface.exchanges,
    interface_biomass = x -> x.interface.biomass,
    default_community_exchange_bound = nothing,
)
    @assert length(models) >= 1 "at least one community member is required"
    @assert growth >= 0 "growth must not be negative"

    bounds_lookup = Dict(community_exchange_bounds)

    constraints = interface_constraints(
        (
            Symbol(k) => if m isa A.AbstractFBCModel
                c = flux_balance_constraints(m; interface)
                (c, interface_exchanges(c), a)
                #TODO scale this
            elseif m isa C.ConstraintTree
                (m, interface_exchanges(m), a)
                # TODO how do we scale this thing most easily?
                # one tree with retained bounds and a few bounds with something else?
                # also how do we scale here with a new variable? :)
            else
                throw(DomainError(m, "unsupported community member type"))
            end for (k, (m, a)) in model_abundances
        );
        out_interface = :community_exchanges,
        out_balance = :community_balance,
        bound = x ->
            get(bounds_lookup, String(last(x)), default_community_exchange_bound),
    )

    growth_sums = [
        Symbol(k) => C.Constraint(sum_value(interface_biomass(constraints[Symbol(k)])))
        for (k, _) in model_abundances
    ]

    constraints *
    :community_balance^C.Constraint(#=TODO sum of Xs=#)
    *
    :community_composition^C.ConstraintTree(
        #=TODO pick out the individual compositions outta the individual trees=#
    )
    
end

export community_composition_balance_constraints

"""
$(TYPEDSIGNATURES)

Run the SteadyCom analysis on several models. All arguments are forwarded to
[`community_composition_balance_constraints`](@ref) which constructs the model;
this function returns the solved model.
"""
community_composition_balance_analysis(args...; kwargs...) = frontend_optimized_values(
    community_composition_balance_constraints,
    args...;
    #TODO zero objective = x -> x.community_biomass.value,
    #TODO feasibility sense
    kwargs...,
)

export community_composition_balance_analysis
