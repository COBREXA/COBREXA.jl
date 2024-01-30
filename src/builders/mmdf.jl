
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

Allocate log-concentration-stoichiometry constraints for the `reaction_subset`
and `metabolite_subset` in `model`, as used by e.g.
[`max_min_driving_force_analysis`](@ref).

The output constraint tree contains a log-concentration variable for each
metabolite in subtree `log_concentrations`. The reaction quotient in log
variables for each reaction is stored in `log_concentration_stoichiometry`.

Function `concentration_bound` may return a bound for the log-concentration of a
given metabolite (compatible with `ConstraintTrees.Bound`), or `nothing`.
"""
function log_concentration_constraints(
    model::A.AbstractFBCModel;
    reaction_subset = A.reactions(model),
    metabolite_subset = A.metabolites(model),
    concentration_bound = _ -> nothing,
)
    stoi = A.stoichiometry(model)

    constraints =
        :log_concentrations^C.variables(
            keys = Symbol.(metabolite_subset),
            bounds = concentration_bound.(metabolite_subset),
        )

    # map old idxs of metabolites to new order of smaller system
    midx_lookup = Dict(
        indexin(metabolite_subset, A.metabolites(model)) .=>
            1:length(metabolite_subset),
    )

    cs = C.ConstraintTree()
    for (ridx_new, ridx_old) in enumerate(indexin(reaction_subset, A.reactions(model)))
        midxs_old, stoich_coeffs = SparseArrays.findnz(stoi[:, ridx_old])
        rid = Symbol(reaction_subset[ridx_new])
        for (midx_old, stoich_coeff) in zip(midxs_old, stoich_coeffs)
            midx_new = midx_lookup[midx_old]
            mid = Symbol(metabolite_subset[midx_new])
            value = constraints.log_concentrations[mid].value * stoich_coeff
            if haskey(cs, rid)
                cs[rid].value += value
            else
                cs[rid] = C.Constraint(; value)
            end
        end
    end

    return constraints * :log_concentration_stoichiometries^cs
end

export log_concentration_constraints
