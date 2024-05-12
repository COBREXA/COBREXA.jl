
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

# # Loading and saving models

using COBREXA

# ## Getting the models reliably from the repositories
#
# For convenience, COBREXA provides a specific function
# [`download_model`](@ref) to download models from repositories that also
# automatically uses the cached downloaded version of the model if it's already
# downloaded, and verifies the checksum to improve reproducibility. It will
# print out a warning in case the model checksum does not match the
# expectation:

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

#md #!!! tip "Tip: How do I get the model hash?"
#md #    We do not need to fill in the hash values immediately -- instead, it's possible to simply run the function once, and then copy the reported hash value from the warning message into the script.

# ## Loading models

# To load genome-scale metabolic models, COBREXA uses the
# [AbstractFBCModels](https://github.com/COBREXA/AbstractFBCModels.jl) framework
# to import various kinds of models including SBML, JSON and the legacy
# Matlab-formatted "COBRA toolbox" models.

# All models can be loaded automatically using [`load_model`](@ref); but one
# must import the model-type specific packages to load the functionality. (This
# step is required to keep the "base" COBREXA as efficient and fast-loading as
# possible.)

import JSONFBCModels, SBMLFBCModels

model1 = load_model("e_coli_core.json")

model2 = load_model("e_coli_core.xml")

# We can now explore the contents of the models using the AbstractFBCModels'
# interface:
import AbstractFBCModels as A

A.reactions(model1)

A.reactions(model2)

# Additional extractable information can be found in [the documentation of the
# abstract models
# package](https://cobrexa.github.io/AbstractFBCModels.jl/stable/reference/#Model-content-accessors).

# ### Converting model types

# Normally, [`load_model`](@ref) is forced to guess the model type from the
# filename suffix. We can specify the model type ourselves (this also allows
# the users to work with non-standard file suffixes):

model = load_model(JSONFBCModels.JSONFBCModel, "e_coli_core.json")

# Sometimes it is useful to convert the model data to another type, such as the
# SBML to a JSON model structure:

model_converted_to_json = load_model("e_coli_core.xml", JSONFBCModels.JSONFBCModel)

# Or to the "Canonical Julia model" from AbstractFBCModels:
model_in_julia_structures =
    load_model(JSONFBCModels.JSONFBCModel, "e_coli_core.json", A.CanonicalModel.Model)

#md #!!! tip "Compatibility with COBREXA v1.x"
#md #    `CanonicalModel` is a newer, cleaned-up version of the `StandardModel` type used in COBREXA version 1. For all code that relied on `StandardModel`, the canonical one should work just as well.

# The above command specifies all model types explicitly, leaving least room
# for guessing-based errors.
#
# If required, it is also possible to convert all model types to each other
# simply by using Julia's `convert`:

model_in_json_structure = convert(JSONFBCModels.JSONFBCModel, model_in_julia_structures)

# ## Saving models

# The models can be saved to file storage by using [`save_model`](@ref):
save_model(model_converted_to_json, "e_coli_core_from_sbml.json")

# Expectably, the file will contain the JSON with the model description:
println(open("e_coli_core_from_sbml.json") do f
    read(f, 100)
end |> String, "...")
