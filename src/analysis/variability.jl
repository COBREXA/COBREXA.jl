
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

Simplified variant of [`constraints_variability`](@ref) that computes the
variability of all values in tree `targets`, and returns a new tree of the same
shape as `targets` that contains tuples for minima and maxima.

All other arguments are forwarded to the matrix-returning overload of
[`constraints_variability`](@ref).
"""
function constraints_variability(
    constraints::C.ConstraintTree,
    targets::C.ConstraintTree;
    kwargs...,
)

    result_array = constraints_variability(
        constraints,
        tree_deflate(C.value, targets, C.Value);
        kwargs...,
    )

    tree_reinflate(
        targets,
        Tuple{eltype(result_array),eltype(result_array)}[
            tuple(a, b) for (a, b) in eachrow(result_array)
        ],
    )
end

"""
$(TYPEDSIGNATURES)

In a feasible space specified by `constraints`, compute the feasible range of
individual `targets` values. The output is a matrix with one column for minima
and second column for maxima of the individual target's values.

This is used e.g. to compute the [`flux_variability_analysis`](@ref), and can
be viewed as a more generalized version thereof.

`output` and `output_type` can be used to customize the information reported
from the solved models.

Extra arguments are passed to [`screen_optimization_model`](@ref).
"""
function constraints_variability(
    constraints::C.ConstraintTree,
    targets::Vector{<:C.Value};
    output = (dir, om) -> dir * J.objective_value(om),
    output_type::Type{T} = Float64,
    kwargs...,
)::Matrix{Maybe{T}} where {T}

    target_array = [(dir, tgt) for tgt in targets, dir in (-1, 1)]

    screen_optimization_model(
        constraints,
        target_array;
        kwargs...
    ) do om, (dir, tgt)
        J.@objective(om, Maximal, C.substitute(dir * tgt, om[:x]))
        J.optimize!(om)
        is_solved(om) ? output(dir, om) : nothing
    end
end

export constraints_variability
