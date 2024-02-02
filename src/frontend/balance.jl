
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

Compute an optimal objective-optimizing solution of the given `model`.

Most arguments are forwarded to [`optimized_values`](@ref).

Returns a tree with the optimization solution of the same shape as
given by [`flux_balance_constraints`](@ref).
"""
flux_balance_analysis(model::A.AbstractFBCModel; kwargs...) = frontend_optimized_values(
    flux_balance_constraints,
    model;
    objective = x -> x.objective.value,
    kwargs...,
)

export flux_balance_analysis
