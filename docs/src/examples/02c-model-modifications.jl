
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

# # Making adjustments to the model
#
# Typically, we do not need to solve the models as they come from the authors
# (someone else already did that!), but we want to perform various
# perturbations in the model structure and conditions, and explore how the
# model behaves in the changed conditions.
#
# With COBREXA, there are 2 different approaches that one can take:
# 1. We can change the model structure, and use the changed metabolic model.
#    This is better for doing simple and small, but systematic modifications,
#    such as removing metabolites, adding reactions, etc.
# 2. We can intercept the pipeline that converts the metabolic model to
#    constraints and to the optimizer representation, and make modifications
#    along that way. This is better suited to making global model adjustments,
#    such as using combined objectives, adding reaction-coupling constraints,
#    and combining multiple models into a bigger one.
#
# Here we demonstrate the first, "modeling" approach. The main advantage of
# this approach is that the modified model is still a FBC model, and we can
# export, save and share it via the AbstractFBCModels interface. The main
# disadvantage is that the "common" FBC model interface does not easily express
# various complicated constructions (communities, reaction coupling, enzyme
# constraints, etc.) -- see the [example about modifying the
# constraints](02d-constraint-modifications.md) for more details.
#
# ## Getting the base model

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels

# For applying the modifications, we will use the canonical model as exported
# from package `AbstractFBCModels`. There are other possibilities, but the
# canonical one is easiest to use for common purposes.

import AbstractFBCModels.CanonicalModel as CM

# We can now load the model:

model = convert(CM.Model, load_model("e_coli_core.json"))

# The canonical model is quite easy to work with, made basically of the most
# accessible Julia structures possible. For example, we can look at a reaction
# as such:

model.reactions["PFK"]

#

model.reactions["CS"].stoichiometry

#md # !!! tip "Create custom model types!"
#md #     For some applications, `CanonicalModel` might be too restrictive. Creating a custom model type that perfectly fits a use-case can be done simply by overloading several functions. The documentation of [AbstractFBCModels](https://github.com/COBREXA/AbstractFBCModels.jl) describes the process closer. Further, all model types that adhere to the AbstractFBCModels' interface will "just work" with _all_ analysis functions in COBREXA!

# ## Running FBA on modified models
#
# Since the canonical model is completely mutable, we can change it in any way
# we like and feed the result directly into [`flux_balance_analysis`](@ref).
# Let's first find a "original" solution, so that we have a base solution for
# comparing:

import HiGHS

base_solution = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
base_solution.objective

# Now, for example, we can limit the intake of glucose by the model:

model.reactions["EX_glc__D_e"]

# Since the original intake limit is 10 units, let's try limiting that to 5:

model.reactions["EX_glc__D_e"].lower_bound = -5.0

# ...and solve the modified model:
#
low_glucose_solution = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)
low_glucose_solution.objective

@test isapprox(low_glucose_solution.objective, 0.41559777, atol = TEST_TOLERANCE) #src

# ## Preventing reference-based sharing problems with `deepcopy`
#
# People often want to try different perturbations with a single base model. It
# would therefore look feasible to save the "unmodified" model in a single
# variable, and make copies of that with the modifications applied. Let's
# observe what happens:

base_model = convert(CM.Model, load_model("e_coli_core.json")) # load the base

modified_model = base_model # seemingly make a "copy" for modification

modified_model.reactions["EX_glc__D_e"].lower_bound = -123.0 # modify the glucose intake limit

# Surprisingly, the base model got modified too!

base_model.reactions["EX_glc__D_e"]

# This is because Julia uses reference-based sharing whenever anything mutable
# is copied using the `=` operator. While this is extremely useful in many
# scenarios for data processing efficiency and computational speed, it
# unfortunately breaks this simple use-case.
#
# To fix this situation, we must always remember to make an actual copy of the
# model data, by either carefully copying the changed parts (e.g., using a
# similar approach as with the "shallow" `copy()`), or simply by copying the
# whole model structure as is with `deepcopy()`. Let's try again:

base_model = convert(CM.Model, load_model("e_coli_core.json"))
modified_model = deepcopy(base_model) # this forces an actual copy of the data
modified_model.reactions["EX_glc__D_e"].lower_bound = -123.0

# With `deepcopy`, the result works as intended:

(
    modified_model.reactions["EX_glc__D_e"].lower_bound,
    base_model.reactions["EX_glc__D_e"].lower_bound,
)

@test modified_model.reactions["EX_glc__D_e"].lower_bound != #src
      base_model.reactions["EX_glc__D_e"].lower_bound #src

#md # !!! danger "Avoid overwriting base models when using in-place modifications"
#md #     Whenever changing a copy of the model, check that the base model is not inadvertently changed via a reference. Always use some copy mechanism such as `copy` or `deepcopy` to prevent the default reference-based sharing.

# ## Observing the differences
#
# We already have a `base_solution` and `low_glucose_solution` from above. What
# is the easiest way to see what has changed? We can quite easily compute
# squared distance between all dictionary entries using Julia function for
# merging dictionaries (called `mergewith`).

# With that, we can extract the plain difference in fluxes:
flux_differences = mergewith(-, base_solution.fluxes, low_glucose_solution.fluxes)

# ...and see what were the biggest directional differences:
sort(collect(flux_differences), by = last)

# ...or compute the squared distance, to see the "absolute" changes:
flux_changes =
    mergewith((x, y) -> (x - y)^2, base_solution.fluxes, low_glucose_solution.fluxes)

# ...and again see what changed most:
sort(collect(flux_changes), by = last)

#md # !!! tip "Always use a uniquely defined flux solutions for flux comparisons"
#md #     Since the usual flux balance allows a lot of freedom in the "solved" flux and the only value that is "reproducible" by the analysis is the objective, one should never compare the flux distributions directly. Typically, that may result in false-positive (and sometimes false-negative) differences. Use e.g. [parsimonious FBA](03b-parsimonious-flux-balance.md) to obtain uniquely determined and safely comparable flux solutions.

# ## Coupling constraints
#
# Some model types support additional constraints over the reaction fluxes,
# which are historically called "coupling". These allow to e.g. place a bound
# on a total flux through several reactions.
#
# Canonical model supports these as "couplings":

model.couplings["total_energy_intake"] = CM.Coupling(
    lower_bound = 0,
    upper_bound = 5,
    reaction_weights = Dict("EX_glc__D_e" => -1.0, "EX_fru_e" => -1.0, "EX_pyr_e" => -1.0),
)

# The values of any coupling constraints can be inspected directly in the
# solved model:

solution_with_coupling = flux_balance_analysis(model, optimizer = HiGHS.Optimizer)

solution_with_coupling.coupling.total_energy_intake

@test 0 <= solution_with_coupling.coupling.total_energy_intake <= 5 #src
@test isapprox(solution_with_coupling.objective, 0.4155977750928965, atol = TEST_TOLERANCE) #src
