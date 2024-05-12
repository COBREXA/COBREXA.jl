
# Copyright (c) 2021-2024, University of Luxembourg                         #src
# Copyright (c) 2021-2024, Heinrich-Heine University Duesseldorf            #src
#                                                                           #src
# Licensed under the Apache License, Version 2.0 (the "License");           #src
# you may not use this file except in compliance with the License.          #src
# You may obtain a copy of the License at                                   #src
#                                                                           #src
#     http://www.apache.org/licenses/LICENSE-2.0                            #src
#                                                                           #src
# Unless required by applicable law or agreed to in writing, software       #src
# distributed under the License is distributed on an "AS IS" BASIS,         #src
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  #src
# See the License for the specific language governing permissions and       #src
# limitations under the License.                                            #src

# # Production envelopes
# Production envelopes determine the flux of the model objective at different
# values of specific reactions, spanning their variability. We can use the
# builtin function [`objective_production_envelope`](@ref) to quickly simulate
# this.

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels
import GLPK

model = load_model("e_coli_core.json")

envelope = objective_production_envelope(
    model,
    ["EX_o2_e", "EX_co2_e"];
    breaks = 5,
    optimizer = GLPK.Optimizer,
)

@test count(isnothing, envelope.objective_values) == 18 #src
