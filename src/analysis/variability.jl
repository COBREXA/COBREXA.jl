
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
function constraints_variability(
    constraints::C.ConstraintTree,
    targets::C.ConstraintTree;
    optimizer,
    settings = [],
    workers = D.workers(),
)::C.Tree{Tuple{Maybe{Float64},Maybe{Float64}}}

    result_array = constraints_variability(
        constraints,
        tree_deflate(C.value, targets, C.Value);
        optimizer,
        settings,
        workers,
    )

    tree_reinflate(
        targets,
        Tuple{Maybe{Float64},Maybe{Float64}}[
            tuple(a, b) for (a, b) in eachrow(result_array)
        ],
    )
end

"""
$(TYPEDSIGNATURES)

TODO
"""
function constraints_variability(
    constraints::C.ConstraintTree,
    targets::Vector{<:C.Value};
    optimizer,
    settings = [],
    output = (dir, om) -> dir * J.objective_value(om),
    output_type::Type{T} = Float64,
    workers = D.workers(),
)::Matrix{Maybe{T}} where {T}

    target_array = [(dir, tgt) for tgt in targets, dir in (-1, 1)]

    screen_optimization_model(
        constraints,
        target_array;
        optimizer,
        settings,
        workers,
    ) do om, (dir, tgt)
        J.@objective(om, Maximal, C.substitute(dir * tgt, om[:x]))
        J.optimize!(om)
        is_solved(om) ? output(dir, om) : nothing
    end
end

export constraints_variability
