
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

# # Making adjustments to the constraint system
#
# In the [previous example about model
# adjustments](02b-model-modifications.md), we noted that some constraint
# systems may be too complex to be changed within the limits of the usual FBC
# model view, and we may require a sharper tool to do the changes we need. This
# example shows how to do that by modifying the constraint systems that are
# generated within COBREXA to represent the metabolic model contents.

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels
import GLPK

model = load_model("e_coli_core.json") # flux balance type model

# ## Background: Constraint trees

# COBREXA uses [ConstraintTrees](https://github.com/COBREXA/ConstraintTrees.jl)
# to represent model structures internally. This framework is incredibly powerful,
# as it neatly groups relevant variables and constraints together.

import ConstraintTrees as C

# In general, constraint-based models use fluxes as variables, and all the
# constraints are in terms of them (or derived quantities). For "normal" models,
# you can directly convert their flux balance format (from json, sbml, mat
# files) into a ConstraintTree. These structures make it particularly easy to
# formulate new constraints.

ctmodel = flux_balance_constraints(model) # load the ConstraintTree of model

# Notice, variables and constraints are grouped here.

ctmodel.fluxes # all variables

#

ctmodel.flux_stoichiometry # mass balance constraints

#

ctmodel.objective # objective (usually specified in the model as a biomass function), notice it does not have a bound

# ## Customizing the model

# ConstraintTrees make is simple to modify the model. The most explicit way to
# do this is to make a new constraint tree representation of the model. But
# first, let's make a new constraint to represent fermentation fluxes.

fermentation = ctmodel.fluxes.EX_ac_e.value + ctmodel.fluxes.EX_etoh_e.value # acetate and ethanol fluxes are grouped

fermentation_constraint = C.Constraint(fermentation, (10.0, 1000.0)) # create a new constraint, bounding the flux

fermentation_constrainttree = :fermentation^fermentation_constraint # create a new ConstraintTree, naming this constraint

forced_mixed_fermentation = ctmodel * fermentation_constrainttree # new modified model is created

# ConstraintTrees can be directly solved. The variables and constraints are
# automatically parsed into a JuMP model, which is subsequently solved. Note,
# you need to specify the objective.

vt = optimized_values(
    forced_mixed_fermentation,
    objective = forced_mixed_fermentation.objective.value,
    optimizer = GLPK.Optimizer,
)

@test isapprox(vt.objective, 0.6337, atol = TEST_TOLERANCE) #src

# Models that cannot be solved return `nothing`. In the example below, the
# underlying model is modified.

ctmodel.fluxes.ATPM.bound = C.Between(1000.0, 10000.0)

vt = optimized_values(
    ctmodel,
    objective = ctmodel.objective.value,
    optimizer = GLPK.Optimizer,
)

@test isnothing(vt) #src

# In general, every attribute of a ConstraintTree can be modified. Using some
# building block functions, complicated models can be formulated. Here we will
# create a model with only positive fluxes, by splitting all the reactions into
# two components (forward and reverse components). This is frequently a first
# step in building more complicated models.

positive_model = deepcopy(ctmodel)

positive_model += sign_split_variables( # notice the +
    positive_model.fluxes,
    positive = :fluxes_forward,
    negative = :fluxes_reverse,
)

#md # !!! warning "Warning: Take care between + and * operators for ConstraintTrees"
#md #    For ConstraintTrees, `+` adds new variables to a model, and `*` adds constraints to a model with the assumption that the variables they reference are already in the model.

# After creating the new variables, we need to link them to the original
# variables, using constraints.

positive_model *=
    :pos_neg_flux_link^sign_split_constraints(; # notice the *
        positive = positive_model.fluxes_forward,
        negative = positive_model.fluxes_reverse,
        signed = positive_model.fluxes,
    )

# Next, we can specify a new objective, minimizing the sum of all positive fluxes.

positive_model *=
    :l1_objective^C.Constraint(
        sum(C.value(v) for v in values(positive_model.fluxes_forward)) +
        sum(C.value(v) for v in values(positive_model.fluxes_reverse)),
        nothing, # no bound
    )

# Notice how easy it was to sum up all the fluxes in the forward and reverse
# directions. Also, we did not lose any information, as the new variables and
# objective are just layered on top of the original model.

# Next, we specify a specific growth rate, as a new constraint (making its
# removal simple later).

positive_model *=
    :growth_rate_setpoint^C.Constraint(
        C.value(positive_model.fluxes.BIOMASS_Ecoli_core_w_GAM),
        C.EqualTo(0.6), # 1/h
    )

l1_sol = optimized_values(
    positive_model,
    objective = positive_model.l1_objective.value,
    optimizer = GLPK.Optimizer,
    sense = COBREXA.Minimal,
)

# Removing constraints is simple.

delete!(positive_model, :l1_objective)

#md # !!! warning "Warning: Take care to keep your model consistent"
#md #    While ConstraintTrees gives you the power to very simply create complex models, it does not guard you against making the internal structure inconsistent (e.g. changing the bounds of the positive variables to allow negative numbers, messing with the link constraints, etc.).
