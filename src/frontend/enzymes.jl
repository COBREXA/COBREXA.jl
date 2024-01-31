
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

"""
$(TYPEDEF)

A simple struct storing information about the isozyme composition, including
subunit stoichiometry and turnover numbers. Use with
[`enzyme_constrained_flux_balance_analysis`](@ref).

# Fields
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct Isozyme
    gene_product_stoichiometry::Dict{String,Float64}
    kcat_forward::Maybe{Float64} = nothing
    kcat_reverse::Maybe{Float64} = nothing
end

export Isozyme

"""
$(TYPEDSIGNATURES)

Construct a enzyme-constrained flux-balance constraint system. The model is
parameterized by `reaction_isozymes`, which is a mapping of reaction
identifiers to [`Isozyme`](@ref) descriptions.

Additionally, the computation requires `gene_product_molar_masses` to describe
the weights of enzyme building material, and `capacity`, which limits the mass
of enzymes in the whole model.

`capacity` may be a single number, which sets the limit for "all described
enzymes". Alternatively, `capacity` may be a vector of identifier-genes-limit
triples that make a constraint (identified by the given identifier) that limits
the listed genes to the given limit.
"""
function enzyme_constrained_flux_balance_constraints(
    model::A.AbstractFBCModel;
    reaction_isozymes::Dict{String,Dict{String,Isozyme}},
    gene_product_molar_masses::Dict{String,Float64},
    capacity::Union{Vector{Tuple{String,Vector{String},Float64}},Float64},
)
    constraints = flux_balance_constraints(model)

    # might be nice to omit some conditionally (e.g. slash the direction if one
    # kcat is nothing)
    isozyme_amounts = isozyme_amount_variables(
        Symbol.(keys(reaction_isozymes)),
        rid -> Symbol.(keys(reaction_isozymes[string(rid)])),
    )

    # allocate variables for everything (nb. += wouldn't associate right here)
    constraints =
        constraints +
        :fluxes_forward^unsigned_positive_contribution_variables(constraints.fluxes) +
        :fluxes_reverse^unsigned_negative_contribution_variables(constraints.fluxes) +
        :isozyme_forward_amounts^isozyme_amounts +
        :isozyme_reverse_amounts^isozyme_amounts +
        :gene_product_amounts^C.variables(
            keys = Symbol.(A.genes(model)),
            bounds = C.Between(0, Inf),
        )

    # connect all parts with constraints
    constraints *
    :directional_flux_balance^sign_split_constraints(
        positive = constraints.fluxes_forward,
        negative = constraints.fluxes_reverse,
        signed = constraints.fluxes,
    ) *
    :isozyme_flux_forward_balance^isozyme_flux_constraints(
        constraints.isozyme_forward_amounts,
        constraints.fluxes_forward,
        (rid, isozyme) -> maybemap(
            x -> x.kcat_forward,
            maybeget(reaction_isozymes, string(rid), string(isozyme)),
        ),
    ) *
    :isozyme_flux_reverse_balance^isozyme_flux_constraints(
        constraints.isozyme_reverse_amounts,
        constraints.fluxes_reverse,
        (rid, isozyme) -> maybemap(
            x -> x.kcat_reverse,
            maybeget(reaction_isozymes, string(rid), string(isozyme)),
        ),
    ) *
    :gene_product_isozyme_balance^gene_product_isozyme_constraints(
        constraints.gene_product_amounts,
        (constraints.isozyme_forward_amounts, constraints.isozyme_reverse_amounts),
        (rid, isozyme) -> maybemap(
            x -> [(Symbol(k), v) for (k, v) in x.gene_product_stoichiometry],
            maybeget(reaction_isozymes, string(rid), string(isozyme)),
        ),
    ) *
    :gene_product_capacity^(
        capacity isa Float64 ?
        C.Constraint(
            value = sum(
                gpa.value * gene_product_molar_masses[String(gp)] for
                (gp, gpa) in constraints.gene_product_amounts
            ),
            bound = C.Between(0, capacity),
        ) :
        C.ConstraintTree(
            Symbol(id) => C.Constraint(
                value = sum(
                    constraints.gene_product_amounts[Symbol(gp)].value *
                    gene_product_molar_masses[gp] for gp in gps
                ),
                bound = C.Between(0, limit),
            ) for (id, gps, limit) in capacity_limits
        )
    )
end

export enzyme_constrained_flux_balance_constraints

"""
$(TYPEDSIGNATURES)

Perform the enzyme-constrained flux balance analysis on the `model` and return the solved constraint system.

Arguments are forwarded to
[`enzyme_constrained_flux_balance_constraints`](@ref); solver configuration
arguments are forwarded to [`optimized_values`](@ref).
"""
enzyme_constrained_flux_balance_analysis(model::A.AbstractFBCModel; kwargs...) =
    frontend_optimized_values(
        enzyme_constrained_flux_balance_constraints,
        model;
        objective = x -> x.objective.value,
        kwargs...,
    )

export enzyme_constrained_flux_balance_analysis
