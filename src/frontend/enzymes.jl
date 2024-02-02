
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

Construct a enzyme-constrained flux-balance constraint system. Based on the
algorithm used in GECKO in *Sánchez, Benjamín J., et al. "Improving the
phenotype predictions of a yeast genome‐scale metabolic model by incorporating
enzymatic constraints." Molecular systems biology 13.8 (2017): 935*.

The model is parameterized by `reaction_isozymes`, which is a mapping of
reaction identifiers to [`Isozyme`](@ref) descriptions. Additionally, the
computation requires `gene_product_molar_masses` to describe the weights of
enzyme building material, and `capacity`, which limits the mass of enzymes in
the whole model.

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
    # these functions should not fail gracefully to prevent user footguns
    # they take strings and but cts use symbols internally, hence lots of conversions are required - kind of annoying
    isozyme_ids(rid) =
        haskey(reaction_isozymes, String(rid)) ?
        Symbol.(keys(reaction_isozymes[String(rid)])) : nothing
    kcat_forward(rid, iso_id) = reaction_isozymes[String(rid)][String(iso_id)].kcat_forward
    kcat_reverse(rid, iso_id) = reaction_isozymes[String(rid)][String(iso_id)].kcat_reverse
    isozyme_subunit_stoichiometry(rid, iso_id) = Dict(
        Symbol(k) => v for (k, v) in
        reaction_isozymes[String(rid)][String(iso_id)].gene_product_stoichiometry
    )
    gene_product_molar_mass(gid) = gene_product_molar_masses[String(gid)]

    capacity_limits =
        capacity isa Real ? [("totalcapacity", Symbol.(A.genes(model)), capacity)] :
        capacity

    enzyme_constraints(
        flux_balance_constraints(model);
        gene_ids = Symbol.(A.genes(model)),
        isozyme_ids,
        kcat_forward,
        kcat_reverse,
        isozyme_subunit_stoichiometry,
        gene_product_molar_mass,
        capacity_limits,
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
