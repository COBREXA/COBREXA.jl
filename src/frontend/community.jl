
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

function community_flux_balance_analysis(
    model_abundances::Vector{Tuple{A.AbstractFBCModel,Float64}},
    optimizer;
    kwargs...,
)
    # TODO f this gets complicated, make a specialized community_constraints
    # builder or so. But ideally this is just module loading + 1 big join +
    # optimizer run.
end
