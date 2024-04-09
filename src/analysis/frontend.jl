
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

A helper that converts a front-end constraint `builder` function (the output of
which would normally be just passed through [`optimized_values`](@ref)) to
front-end analysis function.
"""
function frontend_optimized_values(
    builder,
    args...;
    builder_kwargs = NamedTuple(),
    objective,
    output = identity,
    sense = Maximal,
    optimizer,
    settings = [],
    kwargs...,
)
    constraints = builder(args...; builder_kwargs..., kwargs...)

    # arguments need to be kept in sync
    optimized_values(
        constraints;
        objective = objective(constraints),
        output = output(constraints),
        sense,
        optimizer,
        settings,
    )
end

"""
$(TYPEDSIGNATURES)

A helper that converts a parsimonious-style front-end constraint `builder`
function to front-end analysis function.

Like [`frontend_optimized_values`](@ref), but internally calls
[`parsimonious_optimized_values`](@ref).
"""
function frontend_parsimonious_optimized_values(
    builder,
    args...;
    builder_kwargs = NamedTuple(),
    objective = identity,
    output = identity,
    sense = Maximal,
    optimizer,
    settings = [],
    parsimonious_objective,
    parsimonious_optimizer = nothing,
    parsimonious_sense = Minimal,
    parsimonious_settings = [],
    tolerances = [absolute_tolerance_bound(0)],
    kwargs...,
)
    constraints = builder(args...; builder_kwargs..., kwargs...)

    # arguments need to be kept in sync
    parsimonious_optimized_values(
        constraints;
        objective = objective(constraints),
        output = output(constraints),
        sense,
        settings,
        optimizer,
        parsimonious_objective = parsimonious_objective(constraints),
        parsimonious_optimizer,
        parsimonious_sense,
        parsimonious_settings,
        tolerances,
    )
end
