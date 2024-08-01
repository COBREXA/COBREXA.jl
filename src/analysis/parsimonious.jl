
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

Optimize the system of `constraints` to get the optimal `objective` value. Then
try to find a "parsimonious" solution with the same `objective` value, which
optimizes the `parsimonious_objective` (possibly also switching optimization
sense, optimizer, and adding more settings).

For efficiency, everything is performed on a single instance of JuMP model.

A simpler version suitable for direct work with metabolic models is available
in [`parsimonious_flux_balance_analysis`](@ref).
"""
function parsimonious_optimized_values(
    constraints::C.ConstraintTreeElem;
    objective::C.Value,
    objective_value::Maybe{Float64} = nothing,
    settings = [],
    parsimonious_objective::C.Value,
    parsimonious_optimizer = nothing,
    parsimonious_sense = Minimal,
    parsimonious_settings = [],
    tolerances = [absolute_tolerance_bound(0)],
    output = constraints,
    kwargs...,
)
    # arguments need to be kept in sync with
    # frontend_parsimonious_optimized_values

    om = optimization_model(constraints; objective, kwargs...)
    for m in [configuration.default_solver_settings; settings]
        m(om)
    end

    # first solve the optimization problem with the original objective
    # (if required)
    target_objective_value = if isnothing(objective_value)
        J.optimize!(om)
        is_solved(om) || return nothing
        J.objective_value(om)
    else
        float(objective_value) # make sure it's a pure floaty type
    end

    # switch to parsimonizing the solution w.r.t. to the objective value
    isnothing(parsimonious_optimizer) || J.set_optimizer(om, parsimonious_optimizer)
    for m in [configuration.default_solver_settings; parsimonious_settings]
        m(om)
    end

    J.@objective(om, J.MIN_SENSE, C.substitute(parsimonious_objective, om[:x]))

    # try all admissible tolerances
    for tolerance in tolerances
        (lb, ub) = tolerance(target_objective_value)
        J.@constraint(
            om,
            pfba_tolerance_constraint,
            lb <= C.substitute(objective, om[:x]) <= ub
        )

        J.optimize!(om)
        is_solved(om) && return C.substitute_values(output, J.value.(om[:x]))

        J.delete(om, pfba_tolerance_constraint)
        J.unregister(om, :pfba_tolerance_constraint)
    end

    # all tolerances failed
    return nothing
end

export parsimonious_optimized_values
