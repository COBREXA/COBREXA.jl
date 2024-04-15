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

# # Community FBA models

using COBREXA

# Here we will construct a community FBA model of two  *E. coli* "core" models
# that can interact by exchanging selected metabolites. To do this, we will need
# the model, which we can download if it is not already present.

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

# Additionally to COBREXA and the model format package, we will need a solver
# -- let's use GLPK here:

import JSONFBCModels
import GLPK
import AbstractFBCModels.CanonicalModel as CM

ecoli1 = load_model("e_coli_core.json", CM.Model)
ecoli1.reactions["EX_glc__D_e"].lower_bound = -1000.0
ecoli1.reactions["EX_glc__D_e"].upper_bound = 1000.0
ecoli2 = deepcopy(ecoli1)

# customize models a bit

ecoli1.reactions["CYTBD"].lower_bound = ecoli1.reactions["CYTBD"].upper_bound = 0.0
ecoli2.reactions["FBA"].lower_bound = ecoli2.reactions["FBA"].upper_bound = 0.0

# solve

solution = community_flux_balance_analysis(
    ["bug1" => (ecoli1, 0.2), "bug2" => (ecoli2, 0.8)],
    ["EX_glc__D_e" => (-10.0, 0.0)],
    optimizer = GLPK.Optimizer,
)

@test isapprox(solution.community_biomass, 0.5237157737585179, atol = TEST_TOLERANCE) #src

@test isapprox( #src
    solution.bug1.fluxes.BIOMASS_Ecoli_core_w_GAM, #src
    solution.bug2.fluxes.BIOMASS_Ecoli_core_w_GAM, #src
    atol = TEST_TOLERANCE, #src
) #src

# ## Compare

# We can now e.g. observe the differences in individual pairs of exchanges:

import ConstraintTrees as C

C.zip(
    tuple,
    solution.bug1.interface.exchanges,
    solution.bug2.interface.exchanges,
    Tuple{Float64,Float64},
)

# ## Which composition is best?

screen(0.0:0.1:1.0) do ratio2
    ratio1 = 1 - ratio2
    res = community_flux_balance_analysis(
        ["bug1" => (ecoli1, ratio1), "bug2" => (ecoli2, ratio2)],
        ["EX_glc__D_e" => (-10.0, 0.0)],
        interface = :sbo, # usually more reproducible
        optimizer = GLPK.Optimizer,
    )
    (ratio1, ratio2) => (isnothing(res) ? nothing : res.community_biomass)
end

# ...seems a lot like `bug1` will eventually disappear.

# ## Note about interfaces

# Some work:
flux_balance_constraints(ecoli1, interface = :sbo).interface

# Some generally don't do well:
flux_balance_constraints(ecoli1, interface = :boundary).interface
