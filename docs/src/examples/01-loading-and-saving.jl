
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

# # Loading and saving models!

using COBREXA

# ## Getting the models reliably from the repositories

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.xml",
    "e_coli_core.xml",
    "b4db506aeed0e434c1f5f1fdd35feda0dfe5d82badcfda0e9d1342335ab31116",
)

# TODO: how to get the hashes: specify a dummy value at once and convert them from a warning
#
# TODO: test this :)

# ## Loading models

import JSONFBCModels, SBMLFBCModels

model1 = load_model("e_coli_core.json")

model2 = load_model("e_coli_core.xml")

import AbstractFBCModels as A

A.reactions(model1)

A.reactions(model2)

# ### Converting model types

# Avoid guessing the model type (works also without the suffix in file name):

model = load_model(JSONFBCModels.JSONFBCModel, "e_coli_core.json")

# Directly convert to a different model type

model_converted_to_json = load_model("e_coli_core.xml", JSONFBCModels.JSONFBCModel)

# Or do all at once, load a precisely typed model and convert it to an easily modifiable representation
model_in_julia_structures =
    load_model(JSONFBCModels.JSONFBCModel, "e_coli_core.json", A.CanonicalModel.Model)

# Also possible to convert everything using Julia's `convert`.

# ## Saving models

save_model(model_converted_to_json, "e_coli_core_from_sbml.json")

println(open("e_coli_core_from_sbml.json") do f
    read(f, 100)
end |> String, "...")

# TODO refer to ABCMT docs for more docs
# TODO carefully refer to matlab models
