
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
$(TYPEDSIGNATURES)

Make a `ConstraintTree` that simulates a gene knockout of `knockout_genes` in
the `model` and disables corresponding `fluxes` accordingly.

Keys of the fluxes must correspond to the reaction identifiers in the `model`.

`knockout_genes` may be any collection that support element tests using `in`.
Since the test is done many times, a `Set` is a preferred contained for longer
lists of genes.

All constraints are equality constraints returned in a single flat
`ConstraintTree`.
"""
gene_knockout_constraints(fluxes::C.ConstraintTree, knockout_genes, model::A.AbstractFBCModel) =
    knockout_constraints(fluxes) do rid
        maybemap(
            !,
            A.reaction_gene_products_available(model, string(rid), !in(knockout_genes)),
            false,
        )
    end

"""
$(TYPEDSIGNATURES)

Convenience overload of [`gene_knockout_constraints`](@ref) for knocking out a
single gene (without the necessity to store the gene identifier in a singleton
container).
"""
gene_knockout_constraints(
    fluxes::C.ConstraintTree,
    knockout_gene::String,
    model::A.AbstractFBCModel,
) = gene_knockout_constraints(fluxes, tuple(knockout_gene), model)

export gene_knockout_constraints

"""
$(TYPEDSIGNATURES)

Compute the objective value of the `model` for all knockouts specified by
`gene_combinations`, which is a vector of gene IDs or tuples of gene IDs that
are knocked out in groups.

Returns a vector in the same order as `gene_combinations`.

Extra arguments (mainly, the `optimizer`) are forwarded to
[`screen_optimization_model`](@ref).
"""
function gene_knockouts(
    model::A.AbstractFBCModel,
    gene_combinations::Vector{<:Union{String,NTuple{N,String} where {N}}} = A.genes(model);
    kwargs...,
)
    rxns = A.reactions(model)
    constraints = flux_balance_constraints(model)

    gene_combinations .=> screen_optimization_model(
        constraints,
        gene_combinations;
        objective = constraints.objective.value,
        sense = Maximal,
        kwargs...,
    ) do om, knockout
        con_refs = [
            J.@constraint(om, C.substitute(c.value, om[:x]) == c.bound.equal_to) for
            c in values(gene_knockout_constraints(constraints.fluxes, knockout, model))
        ]
        res = optimized_objective(om)
        J.delete.(Ref(om), con_refs)
        return res
        # TODO the above construction of con_refs is shared with envelopes and
        # might be also reusable elsewhere. Also also the @constraint creation
        # is kinda similar to what solver.jl does.
        #
        # Makes a nice code saving opportunity.
    end
end

export gene_knockouts
