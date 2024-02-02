
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

Build log-concentration-stoichiometry constraints for the `model`, as used e.g.
by [`max_min_driving_force_analysis`](@ref).

The output constraint tree contains a log-concentration variable for each
metabolite in subtree `log_concentrations`. Individual reactions' total
reactant log concentrations (i.e., all log concentrations of actual reactants
minus all log concentrations of products) have their own variables in
`reactant_log_concentrations`.

Function `concentration_bound` may return a bound for the log-concentration of
a given metabolite (compatible with `ConstraintTrees.Bound`), or `nothing`.
"""
function log_concentration_constraints(
    model::A.AbstractFBCModel;
    concentration_bound = _ -> nothing,
)
    rxns = Symbol.(A.reactions(model))
    mets = Symbol.(A.metabolites(model))
    stoi = A.stoichiometry(model)

    constraints =
        :log_concentrations^C.variables(keys = mets, bounds = concentration_bound.(mets))

    cs = C.ConstraintTree()

    for (midx, ridx, coeff) in zip(SparseArrays.findnz(stoi)...)
        rid = rxns[ridx]
        value = constraints.log_concentrations[mets[midx]].value * coeff
        if haskey(cs, rid)
            cs[rid].value += value
        else
            cs[rid] = C.Constraint(; value)
        end
    end

    return constraints * :reactant_log_concentrations^cs
end

export log_concentration_constraints
