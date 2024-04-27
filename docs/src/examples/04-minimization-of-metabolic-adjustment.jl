
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

# # Minimization of metabolic adjustment analysis

# Minimization of metabolic adjustment analysis (MOMA) finds a flux solution
# that is closest to some reference solution. This may correspond to realistic
# adjustment of living organisms to various perturbations, such as gene
# knockout or environmental stress.
#
# To demonstrate, let's use the E. coli model.

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)


# We shall use both quadratic and linear solvers -- the "closest to some
# reference solution" typically refers to Euclidean ("L2") distance which
# requires a QP solver, but Manhattan ("L1") distance is also demonstrated
# below.
import Clarabel, GLPK

# Because we will have to perform some perturbation, we import the model in
# canonical Julia structures:
import JSONFBCModels
import AbstractFBCModels.CanonicalModel as CM
ecoli = load_model("e_coli_core.json", CM.Model)

# This will be a good reaction for perturbing:
ecoli.reactions["CYTBD"]

# To do the perturbation, we create a model of a strain which has mild issues
# with running the CYTBD reaction. We use `deepcopy` to completely avoid any
# reference sharing issues.
limited_ecoli = deepcopy(ecoli)
limited_ecoli.reactions["CYTBD"].upper_bound = 10.0

# ## Finding parsimonious solutions

# Becuase we are interested in realistic flux distributions, we have to use an
# analysis method which gives one -- in this case, the parsimonious FBA will do
# just right. For later comparison, we first get the optimal parsimonious flux
# distribution in the perturbed model:

solution = parsimonious_flux_balance_analysis(
    limited_ecoli,
    optimizer = Clarabel.Optimizer,
    settings = [silence],
)

# Now, how much is the flux going to differ if we assume the bacterium did only
# minimal adjustment from the previous state with unlimited CYTBD?

moma_solution = metabolic_adjustment_minimization_analysis(
    limited_ecoli, # the model to be examined
    ecoli; # the model that gives the reference flux
    optimizer = Clarabel.Optimizer,
    settings = [silence],
)

@test isapprox(moma_solution.objective, 0.241496699165187, atol = TEST_TOLERANCE) #src
@test sqrt(moma_solution.minimal_adjustment_objective) < 71 #src
@test isapprox(moma_solution.fluxes.CYTBD, 10.0, atol = TEST_TOLERANCE) #src

# ## Comparing the results
#
# The difference between the naive and minimally-adjusting solutions can be
# extracted using the constraint tree functionality:

import ConstraintTrees as C
difference = C.zip(-, solution, moma_solution, Float64)

#

sort(collect(difference.fluxes), by = last)

@test isapprox( #src
    maximum(last.(difference.fluxes) .^ 2), #src
    438.56051577031155, #src
    atol = TEST_TOLERANCE, #src
) #src

# ## Using a custom reference flux
#
# In certain situations, you might want to examine how the model would adjust
# from a known reaction flux. You can supply it
# manually as the second argument (instead of the reference model).

ref = parsimonious_flux_balance_analysis(
    ecoli,
    optimizer = Clarabel.Optimizer,
    settings = [silence],
)

ref_closest_solution = metabolic_adjustment_minimization_analysis(
    limited_ecoli,
    ref.fluxes;
    optimizer = Clarabel.Optimizer,
    settings = [silence],
)

@test isapprox( #src
    ref_closest_solution.objective, #src
    moma_solution.objective, #src
    atol = TEST_TOLERANCE, #src
) #src
@test C.reduce( #src
    max, #src
    init = 0, #src
    C.zip((a, b) -> (a - b)^2, ref_closest_solution, moma_solution, Float64), #src
) < TEST_TOLERANCE #src

# The flux may even be partial (which is common with measured fluxes):

measured_fluxes =
    C.Tree{Float64}(:EX_ac_e => 5.0, :EX_o2_e => -2.0, :BIOMASS_Ecoli_core_w_GAM => 0.7)

solution_close_to_measurement = metabolic_adjustment_minimization_analysis(
    limited_ecoli,
    measured_fluxes;
    optimizer = Clarabel.Optimizer,
    settings = [silence],
)

@test isapprox(solution_close_to_measurement.objective, 0.272247, atol = TEST_TOLERANCE) #src

# ## Efficient linear-metric MOMA
#
# The linear version of MOMA avoids having to use the quadratic optimizer in
# the process, giving you more optimizer choices and (typically) much better
# performance. Linear MOMA has the same interface as the quadratic one:

linear_moma_solution = linear_metabolic_adjustment_minimization_analysis(
    limited_ecoli,
    ecoli;
    optimizer = GLPK.Optimizer,
    settings = [silence],
)

sort(collect(linear_moma_solution.fluxes), by = last)

# How much does the flux distribution differ from the L2 solution?

sort(
    collect(C.zip(-, linear_moma_solution.fluxes, moma_solution.fluxes, Float64)),
    by = last,
)
