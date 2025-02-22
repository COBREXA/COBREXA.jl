
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

Perform a binary interval-splitting scan on the interval from `base` to
`extreme`, trying to find a value closest to `extreme` where the constraints
returned by function `constraints_at_value` for a given numerical value between
`base` and `extreme` still contain feasible solutions. The interval splitting
continues until the size of the interval becomes smaller than `tolerance`; at
that point the last solution is returned.

Arguments `optimizer` and `settings` are passed to [`optimized_value`](@ref).
`output` is a function that restricts the value-specific constraint tree to a
set of constraints which should be substitued into and returned in the
resulting value tree.

If no feasible solution is found, the function returs `nothing`.
"""
function feasibility_threshold(
    constraints_at_value,
    base::Float64,
    extreme::Float64;
    tolerance::Float64,
    optimizer,
    settings = [],
    output = identity,
)

    best = nothing

    while abs(base - extreme) >= tolerance
        target = (base + extreme) / 2
        cs = constraints_at_value(target)
        x = optimized_value(
            constraints;
            optimizer,
            settings,
            sense = Feasible,
            output = output(cs),
        )

        if isnothing(x)
            extreme = target
        else
            base = target
            best = x
        end
    end

    return best
end
