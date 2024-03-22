
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

TODO
"""
function gene_knockouts(
    model::A.AbstractFBCModel,
    gene_combinations::Vector{<:Union{String,NTuple{N,String} where {N}}} = A.genes(model);
    workers = D.workers(),
    optimizer,
    settings = [],
    kwargs...,
)
    rxns = A.reactions(model)
    constraints = flux_balance_constraints(model)

    gene_combinations .=> screen_optimization_model(
        constraints,
        gene_combinations;
        objective = constraints.objective.value,
        sense = Maximal,
        optimizer,
        settings,
        workers,
    ) do om, knockout
        con_refs = [
            J.@constraint(om, C.substitute(c.value, om[:x]) == c.bound.equal_to) for
            c in values(knockout_constraints(constraints.fluxes, knockout, model))
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
