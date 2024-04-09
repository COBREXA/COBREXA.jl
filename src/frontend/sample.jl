
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

Run a sampling algorithm on the near-optimal feasible space of the `model` (as
specified by `objective_bound`). By default, the sampling algorithm is ACHR
(the `method` parameter is defaulted to [`sample_chain_achr`](@ref)). The
sampling algorithm is ran for `n_chains` and the iterations for collecting the
sampled values are specified by `collect_iterations`.

`optimizer` is used to generate the warmup (with `settings`) for the sampler
using the usual unidimensional maximum-variability fluxes (as from
[`flux_variability_analysis`](@ref)). All computations are parallelized across
`workers`.

Extra arguments are forwarded to [`sample_constraints`](@ref). Eventually the
arguments will reach the `method` function, so extra arguments can be also used
to customize the methods (e.g., by setting the `epsilon` for the ACHR sampler).
"""
function flux_sample(
    model::A.AbstractFBCModel;
    seed::UInt = rand(UInt),
    objective_bound = relative_tolerance_bound(0.9),
    method = sample_chain_achr,
    n_chains = 10,
    collect_iterations = [32],
    optimizer,
    settings = [],
    workers = D.workers(),
    kwargs...,
)

    constraints = flux_balance_constraints(model)

    objective = constraints.objective.value

    objective_flux = optimized_values(
        constraints;
        objective = constraints.objective.value,
        output = constraints.objective,
        optimizer,
        settings,
    )

    isnothing(objective_flux) && return nothing

    constraints *= :objective_bound^C.Constraint(objective, objective_bound(objective_flux))

    warmup = vcat(
        (
            transpose(v) for (_, vs) in constraints_variability(
                constraints,
                constraints.fluxes;
                optimizer,
                settings,
                output = (_, om) -> J.value.(om[:x]),
                output_type = Vector{Float64},
                workers,
            ) for v in vs
        )...,
    )

    sample_constraints(
        method,
        constraints;
        seed,
        output = constraints.fluxes,
        start_variables = warmup,
        n_chains,
        collect_iterations,
        workers,
        kwargs...,
    )
end

export flux_sample
