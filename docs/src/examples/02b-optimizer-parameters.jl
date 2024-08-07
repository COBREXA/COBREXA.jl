
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

# # Changing optimizer parameters
#
# Many optimizers require fine-tuning to produce best results. We can pass in
# additional optimizer settings via the `settings` parameter of
# [`flux_balance_analysis`](@ref). These include e.g.
#
# - [`set_optimizer_attribute`](@ref), allowing us to tune e.g. iteration
#   limits, tolerances, or floating-point precision (see JuMP documentation for
#   more solver-specific settings)
# - [`set_objective_sense`](@ref), allowing the user to change and reverse the
#   optimization direction, if required
# - [`silence`](@ref) for disabling the debug output of the optimizers
# - [`set_optimizer`](@ref) for replacing the optimizer implementation used
#   (this is not quite useful in this case, but becomes beneficial with more
#   complex, multi-stage optimization problems)
# - [`set_time_limit`](@ref) for putting a time limit on the solver
#   computation (this is quite useful for MILP solvers)
#
# To demonstrate this, let's use the usual toy model:

using COBREXA
import JSONFBCModels, Tulip

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

model = load_model("e_coli_core.json")

# Running a FBA with a silent optimizer that has slightly increased iteration
# limit for IPM algorithm may now look as follows:

solution = flux_balance_analysis(
    model,
    optimizer = Tulip.Optimizer,
    settings = [set_optimizer_attribute("IPM_IterationsLimit", 1000)],
)

@test !isnothing(solution) #src

# To see some of the effects of the configuration changes, we may e.g.
# deliberately cripple the optimizer's possibilities to a few iterations and
# only a little time, which will cause it to fail and return no solution:

solution = flux_balance_analysis(
    model,
    optimizer = Tulip.Optimizer,
    settings = [set_optimizer_attribute("IPM_IterationsLimit", 2), set_time_limit(0.1)],
)

println(solution)

@test isnothing(solution) #src

# To see what failed, users may examine the solver output. Because all solver
# output is silenced by default for efficiency reasons, we need to explicitly
# pass in the [`unsilence`](@ref) setting:

solution = flux_balance_analysis(
    model,
    optimizer = Tulip.Optimizer,
    settings = [
        set_optimizer_attribute("IPM_IterationsLimit", 2),
        set_time_limit(0.1),
        unsilence,
    ],
)

# Applicable optimizer attributes are documented in the documentations of the
# respective optimizers. To browse the possibilities, one might want to see the
# [JuMP documentation page that summarizes the references to the available
# optimizers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers).
#
# Default solver settings can be examined and changed via
# [`Configuration`](@ref).
