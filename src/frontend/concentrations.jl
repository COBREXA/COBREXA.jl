
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

Return log-concentration-stoichiometry constraints for the `model`, as used e.g.
by [`max_min_driving_force_analysis`](@ref).

The output constraint tree contains a log-concentration variable for each
metabolite in `model` if `subset_metabolites` is true for the respective
metabolite id. The substree is called `log_concentrations`. The output
constraint tree also contains a reaction metabolite stoichiometry constraints
for reach reaction in `model` `subset_reactions` is true for the respective
reaction id. The substree is called `reaction_stoichiometry`. 

The keyword arguments are functions:
- `subset_reactions` takes a reaction ID as its argument, and returns true if
  the reaction should be included in the `reaction_stoichiometry` substree.
- `subset_metabolites` takes a metaboilte ID as its argument, and returns true
  if the metabolite should be allocated as a variable in `log_concentrations`.
- `concentration_bound` takes a metabolite ID as its argument, may return either
  a bound for the log-concentration of a given metabolite (compatible with
  `ConstraintTrees.Bound`), or `nothing` if no bound is associated with the
  metabolite log concentration.
"""
function log_concentration_constraints(
    model::A.AbstractFBCModel;
    subset_reactions = _ -> true,
    subset_metabolites = _ -> true,
    concentration_bound = _ -> nothing,
)
    rxns = Dict(
        idx => rid for (idx, rid) in enumerate(A.reactions(model)) if subset_reactions(rid)
    )
    mets = Dict(
        idx => mid for
        (idx, mid) in enumerate(A.metabolites(model)) if subset_metabolites(mid)
    )
    stoi = A.stoichiometry(model)

    # these are the only variables
    constraints =
        :log_concentrations^C.variables(
            keys = Symbol.(values(mets)),
            bounds = concentration_bound.(values(mets)),
        )

    cs = C.ConstraintTree()

    for (midx, ridx, coeff) in zip(SparseArrays.findnz(stoi)...)
        haskey(rxns, ridx) || continue
        rid = Symbol(rxns[ridx])
        mid = Symbol(mets[midx])
        value = constraints.log_concentrations[mid].value * coeff
        if haskey(cs, rid)
            cs[rid].value += value
        else
            cs[rid] = C.Constraint(; value)
        end
    end

    return constraints * :reaction_stoichiometry^cs
end

export log_concentration_constraints
