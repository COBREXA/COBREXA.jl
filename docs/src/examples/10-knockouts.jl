
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

# # Gene knockouts
# FBA is classically very good at predicting the effect of knocking out genes in
# an organism.

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels
import GLPK

model = load_model("e_coli_core.json")

# ## Single gene knockouts

ko_objective_values = gene_knockouts(model, optimizer = GLPK.Optimizer)

ko_dict = Dict(ko_objective_values)

@test length(ko_dict) == 137 #src

ko_dict["b3919"]

ko_dict["b3738"]

# how many genes are critical?

critical = count(isnothing, values(ko_dict))

@test critical == 2 #src

# ## Multiple gene knockouts

some_double_knockouts = gene_knockouts(
    model,
    [("b3919", "b3738"), ("b0118", "b0720")],
    optimizer = GLPK.Optimizer,
)

@test isapprox(last(some_double_knockouts[1]), 0.13475540327383498, atol = TEST_TOLERANCE) #src
@test isapprox(last(some_double_knockouts[2]), 0, atol = TEST_TOLERANCE) #src

knockouts_with_b3919 = gene_knockouts(
    model,
    tuple.(keys(ko_dict), "b3919"),
    optimizer = GLPK.Optimizer,
    settings = [silence],
)

# how many genes are critical with b3919 missing?

critical_without_b3919 = count(isnothing, last.(knockouts_with_b3919))

@test critical_without_b3919 == 23 #src
