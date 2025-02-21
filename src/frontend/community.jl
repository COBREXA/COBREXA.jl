
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
    community_exchange_bounds = Dict();
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

    # TODO: this should be prefixed by `:community` to avoid name collisions.
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

Build a community of `models` of fixed total biomass `growth`, where the
community abundances are represented as variables. This allows finding feasible
community compositions at a given growth rate, following the methodology of
SteadyCom algorithm (Chan SH, Simons MN, Maranas CD. *SteadyCom: predicting
microbial abundances while ensuring community stability*. PLoS computational
biology. 2017 May 15;13(5):e1005539).

`models` is an interable of pairs in the shape of `(name, model)`, where a
`model` can be a constraint tree or an abstract metabolic model. The models are
connected by interfaces and bounded just like in
`community_flux_balance_constraints` (as the main difference, this function
does not accept abundances, but requires knowledge of growth in advance).

The output constraint tree contains the target growth value as `growth`,
relative abundances in subtree `abundances`, all values of the original
models in subtree `community`, and SteadyCom-style diluted constraints in
subtree `diluted_constraints`.
"""
function community_composition_balance_constraints(
    models,
    growth,
    community_exchange_bounds = Dict();
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
                (:community^c, interface_exchanges(c))
            elseif m isa C.ConstraintTree
                (:community^m, interface_exchanges(m))
            else
                throw(DomainError(m, "unsupported community member type"))
            end for (k, m) in models
        );
        out_interface = :community_exchanges,
        out_balance = :community_balance,
        bound = x ->
            get(bounds_lookup, String(last(x)), default_community_exchange_bound),
    )

    constraints +=
        :abundances^C.variables(
            keys = keys(constraints.community),
            bounds = C.Between(0, 1),
        )

    return C.ConstraintTree(
        :community => erase_bounds(constraints.community),
        :diluted_constraints => C.ConstraintTree(
            k => value_scaled_bound_constraints(v, constraints.abundances[k]) for
            (k, v) in constraints.community
        ),
        :abundances => constraints.abundances,
        :total_abundance_constraint => C.Constraint(sum_value(constraints.abundances), 1.0),
        :growth => C.Constraint(C.LinearValue(growth)),
        :biomass_constraints => C.ConstraintTree(
            k => equal_value_constraint(
                sum_value(interface_biomass(constraints.community[k])),
                constraints.abundances[k].value * growth,
            ) for (k, v) in constraints.community
        ),
    )
end

export community_composition_balance_constraints

"""
$(TYPEDSIGNATURES)

Run a SteadyCom-style analysis on a community formed by given models, scanning
for the maximum achievable growth. The scanning proceeds sequentially by
halving the interval `(0, maximum_growth)`, until the interval size becomes
smaller than the desired `tolerance`.

Positional and keyword arguments are forwarded to
[`community_composition_balance_constraints`](@ref). `optimizer` and `settings`
are forwarded to [`optimization_model`](@ref).

Function `output` may be specified to restrict the reported result by reducing
the best-growing constraint tree to a sensible set of constraints to report.
`output` is executed on the tree before the solution variables are substituted
in. For example to just return the best achieved growth and abundances, use
`output = x -> :growth^x.growth * :abundances^x.abundances`.
"""
function community_composition_balance_analysis(
    models,
    maximum_growth,
    args...;
    tolerance = 1e-4,
    optimizer,
    settings = [],
    output = identity,
    kwargs...,
)
    @assert maximum_growth >= 0 "maximum_growth must not be negative"

    (gl, gu) = (0, maximum_growth)

    best_solution = nothing

    while gu - gl >= tolerance
        g = (gu - gl) / 2
        constraints =
            community_composition_balance_constraints(models, g, args...; kwargs...)
        m = optimized_model(
            optimization_model(constraints, sense = Feasible; optimizer, settings),
        )
        if isnothing(m)
            gu = g
        else
            gl = g
            best_solution = (variable_vector(m), constraints)
        end
    end

    if isnothing(best_solution)
        return nothing
    else
        (vars, cs) = best_solution
        return C.substitute_values(output(cs), vars)
    end
end

export community_composition_balance_analysis

"""
$(TYPEDSIGNATURES)

Run a SteadyCom-style analysis, and return the variability of the community
composition at the given growth rate.

Positional and keyword arguments are forwarded to
[`community_composition_balance_constraints`](@ref); `optimizer`, `settings`
and `workers` are internally forwarded to [`constraints_variability`](@ref).
"""
function community_composition_variability_analysis(
    models,
    growth,
    args...;
    optimizer,
    settings = [],
    workers = D.workers(),
    kwargs...,
)
    constraints =
        community_composition_balance_constraints(models, growth, args...; kwargs...)
    return constraints_variability(
        constraints,
        constraints.abundances;
        optimizer,
        workers,
        settings,
    )
end

export community_composition_variability_analysis
