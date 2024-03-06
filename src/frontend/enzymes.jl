
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
    """
    Mapping of gene product identifiers ("genes" in FBC model nomenclature)
    to their relative amount required to construct one unit of the isozyme.
    """
    gene_product_stoichiometry::Dict{String,Float64}

    "Turnover number for this isozyme catalyzing the forward direction of the
    reaction."
    kcat_forward::Maybe{Float64} = nothing

    "Turnover number for this isozyme catalyzing the reverse direction of the
    reaction."
    kcat_reverse::Maybe{Float64} = nothing
end

export Isozyme

"""
$(TYPEDSIGNATURES)

Construct a enzyme-constrained flux-balance constraint system,
following the method in GECKO algorithm (refer to: *Sánchez, Benjamín J., et
al. "Improving the phenotype predictions of a yeast genome‐scale metabolic
model by incorporating enzymatic constraints." Molecular systems biology 13.8
(2017): 935*).

The enzyme mass constraints depend primarily on the available *isozymes*, given
in parameter `reaction_isozymes`, which is a mapping of reaction identifiers to
descriptions of [`Isozyme`](@ref)s that may catalyze the particular reactions.
The isozymes are built from gene products, the mass of which is specified by
`gene_product_molar_masses`. In total, the amount of gene product building
material is limited by `capacity`.

`capacity` may be a single number, which sets the mass limit for "all described
enzymes". Alternatively, `capacity` may be a vector of identifier-genes-limit
triples that together form a constraint (identified by the given identifier)
that limits the total sum of the listed genes to the given limit.

`interface` and `interface_name` are forwarded to
[`flux_balance_constraints`](@ref).
"""
function enzyme_constrained_flux_balance_constraints(
    model::A.AbstractFBCModel;
    reaction_isozymes::Dict{String,Dict{String,Isozyme}},
    gene_product_molar_masses::Dict{String,Float64},
    capacity::Union{Vector{Tuple{String,Vector{String},Float64}},Float64},
    interface::Maybe{Symbol} = nothing,
    interface_name = :interface,
)
    # prepare some accessor functions for the later stuff
    # TODO: might be nicer to somehow parametrize the fwd/rev directions out.
    # Also there is a lot of conversion between symbols and strings, might be
    # nicer to have that sorted out in some better way.
    function isozyme_forward_ids(rid)
        haskey(reaction_isozymes, String(rid)) || return nothing
        return [
            Symbol(k) for
            (k, i) in reaction_isozymes[String(rid)] if !isnothing(i.kcat_forward)
        ]
    end
    function isozyme_reverse_ids(rid)
        haskey(reaction_isozymes, String(rid)) || return nothing
        return [
            Symbol(k) for
            (k, i) in reaction_isozymes[String(rid)] if !isnothing(i.kcat_reverse)
        ]
    end
    kcat_forward(rid, iso_id) = reaction_isozymes[String(rid)][String(iso_id)].kcat_forward
    kcat_reverse(rid, iso_id) = reaction_isozymes[String(rid)][String(iso_id)].kcat_reverse
    isozyme_gene_product_stoichiometry(rid, iso_id) = Dict(
        Symbol(k) => v for (k, v) in
        reaction_isozymes[String(rid)][String(iso_id)].gene_product_stoichiometry
    )
    gene_ids = Symbol.(keys(gene_product_molar_masses))
    gene_product_molar_mass(gid) =
        let k = String(gid)
            if haskey(gene_product_molar_masses, k)
                gene_product_molar_masses[k]
            else
                throw(DomainError(k, "missing a required gene product molar mass"))
                # TODO: maybe better default to 0 and ignore the issue ?
            end
        end

    # allocate all variables and build the system
    constraints = flux_balance_constraints(model; interface, interface_name)

    constraints += sign_split_variables(
        constraints.fluxes,
        positive = :fluxes_forward,
        negative = :fluxes_reverse,
    )

    constraints += enzyme_variables(;
        fluxes_forward = constraints.fluxes_forward,
        fluxes_reverse = constraints.fluxes_reverse,
        gene_ids,
        isozyme_forward_ids,
        isozyme_reverse_ids,
    )

    return constraints *
           sign_split_constraints(;
               positive = constraints.fluxes_forward,
               negative = constraints.fluxes_reverse,
               signed = constraints.fluxes,
           ) *
           enzyme_constraints(;
               fluxes_forward = constraints.fluxes_forward,
               fluxes_reverse = constraints.fluxes_reverse,
               isozyme_forward_amounts = constraints.isozyme_forward_amounts,
               isozyme_reverse_amounts = constraints.isozyme_reverse_amounts,
               gene_product_amounts = constraints.gene_product_amounts,
               kcat_forward,
               kcat_reverse,
               isozyme_gene_product_stoichiometry,
               gene_product_molar_mass,
               capacity_limits = capacity isa Real ?
                                 [(:total_capacity, gene_ids, capacity)] :
                                 [(k, Symbol.(gs), cap) for (k, gs, cap) in capacity],
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
