
# Copyright (c) 2024, University of Luxembourg                              #src
# Copyright (c) 2024, Heinrich-Heine University Duesseldorf                 #src
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

# # Gap filling
#
# Gapfilling is used to find easiest additions to the models that would make
# them feasible and capable of growth.
#
# Typically, an infeasible model ("with gaps") is used together with an
# universal model (which contains "everything"), and the algorithm attempts to
# find the minimal amount of reactions from the universal model that make the
# gappy model happy. In turn, the gapfilling optimization problem becomes a
# MILP.
#
# Gapfilling is sometimes used to produce "viable" genome-scale
# reconstructions from partial ones, but without additional manual intervention
# the gapfilling results are almost never biologically valid. A good use of
# gapfilling is to find problems in a model that cause infeasibility as
# follows: First the modeller makes a set of (unrealistic) universal reactions
# that supply or remove metabolites, and after gapfilling, metabolites that had
# to be supplied or removed to make the model feasible mark possible problems,
# thus guiding the modeller towards correct solution.

# We will use a partially crippled *E. coli* toy model and see the minimal
# amount of reactions that may save it.

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

import JSONFBCModels, HiGHS
model = load_model("e_coli_core.json")

# First, let's produce an infeasible model:

import AbstractFBCModels.CanonicalModel as CM
infeasible_model = convert(CM.Model, model)

for rxn in ["TALA", "PDH", "PGI", "PYK"]
    infeasible_model.reactions[rxn].lower_bound = 0.0
    infeasible_model.reactions[rxn].upper_bound = 0.0
end

# After removing the above reactions, the model will fail to solve:

flux_balance_analysis(infeasible_model, optimizer = HiGHS.Optimizer) |> println

# To avoid very subtle semantic issues, we are going to remove the ATP
# maintenance pseudoreaction from the universal model:
universal_model = convert(CM.Model, model)
delete!(universal_model.reactions, "ATPM")

# ## Making the model feasible with a minimal set of reactions

# Which of the reactions we have to fill back to get the model working again?
# First, let's run [`gap_filling_analysis`](@ref) to get a solution for a
# system that implements the reaction patching:

x = gap_filling_analysis(
    infeasible_model,
    universal_model,
    0.05,
    optimizer = HiGHS.Optimizer,
)

# The reactions that had to be re-added can be found from the `fill_flags`:

filled_reactions = [k for (k, v) in x.fill_flags if v != 0]

@test length(filled_reactions) == 1 #src

# If we want to try to generate another solution, we have to explicitly ask the
# system to avoid the previous solution. That is done via setting the argument
# `known_fill`. We can also set the `max_cost` to avoid finding too benevolent
# solutions:

x2 = gap_filling_analysis(
    infeasible_model,
    universal_model,
    0.05,
    max_cost = 2.0,
    known_fills = [x.fill_flags],
    optimizer = HiGHS.Optimizer,
)

other_filled_reactions = [k for (k, v) in x2.fill_flags if v != 0]

#md # !!! warning "Why is the gapfilling algorithm adding seemingly unneeded reactions?"
#md #     By default, COBREXA does not do any "cleaning" on the universal model; all reactions that are present in that model will be potentially utilized in the new model, and all of them will need to respect their original bounds in the universal model. That becomes an issue with **reactions that are bounded to non-zero flux** (such as the `ATPM` reaction in the E. coli "core" model) -- since their flux is marked as necessarily non-zero to make any model feasible; they will need to be in the fill set, because otherwise their flux would be equal zero.
#md #
#md #     As the simplest solution, all realistic uses of gapfilling should carefully check the set of universal reactions, and ideally exclude all exchanges and pseudoreactions.

# ## Model debugging: which metabolite is missing?
#
# Gap-filling is great for detecting various broken links and imbalances in
# metabolic models. We show how to find the metabolites are causing the
# imbalance for our "broken" E. coli model.
#
# First, we construct a few completely unnatural reactions that create/remove
# the metabolites from/to nowhere:

magic_model = convert(CM.Model, model)
empty!(magic_model.genes)
empty!(magic_model.reactions)

for mid in keys(magic_model.metabolites)
    magic_model.reactions[mid] = CM.Reaction(
        lower_bound = -100.0,
        upper_bound = 100.0,
        stoichiometry = Dict(mid => 1.0),
    )
end

# Gapfilling now points to the metabolites that need to be somehow taken care
# of by the modeller in order for the model to become feasible:

xm = gap_filling_analysis(infeasible_model, magic_model, 0.05, optimizer = HiGHS.Optimizer)

blocking_metabolites = [k for (k, v) in xm.fill_flags if v != 0]

@test length(blocking_metabolites) == 1 #src

# We can also have a look at how much of a given metabolite was used to make
# the model feasible again:

xm.universal_fluxes[first(blocking_metabolites)]
