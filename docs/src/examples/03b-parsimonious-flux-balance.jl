
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

# # Parsimonious flux balance analysis

# Here, we use [`parsimonious_flux_balance_analysis`](@ref) (pFBA) to find the
# optimal flux distribution in the *E. coli* "core" model. In essence, pFBA
# first uses FBA to find an optimal objective value for the model, and then
# minimizes the squared distance of the flux from the zero (i.e., minimizes its
# L2 norm). As the main benefit, this gives a unique (and possibly more
# realistic) solution to the model.

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

# Notably, we need an optimizer that can solve quadratic (QP) models:
import Clarabel

import JSONFBCModels
model = load_model("e_coli_core.json")

# The pFBA is, in its most default form, implemented in function
# [`parsimonious_flux_balance_analysis`](@ref):

solution = parsimonious_flux_balance_analysis(model; optimizer = Clarabel.Optimizer)

#

solution.fluxes

@test isapprox(solution.objective, 0.873922; atol = TEST_TOLERANCE) #src

@test isapprox( #src
    solution.parsimonious_objective, #src
    11414.211988071253, #src
    atol = QP_TEST_TOLERANCE, #src
) #src

@test isapprox( #src
    sum(x^2 for x in values(solution.fluxes)), #src
    solution.parsimonious_objective, #src
    atol = QP_TEST_TOLERANCE, #src
) #src

# ## Using different solvers for the problem stages
#
# Sometimes it is useful to employ a dedicated LP solver to find the solution to
# the original FBA problem, and then a dedicated QP solver to minimize the
# fluxes. We can set the optimizer and parsimonious optimizer separately using
# keyword arguments:

import HiGHS

solution = parsimonious_flux_balance_analysis(
    model;
    optimizer = HiGHS.Optimizer, # HiGHS is used only for LP here
    parsimonious_optimizer = Clarabel.Optimizer, # Clarabel is great for solving QPs
)

@test isapprox(solution.objective, 0.873922; atol = TEST_TOLERANCE) #src

# ## Using linear (L1) parsimonious constraints

# For efficiency reasons, it is also possible to use a pFBA version that
# optimizes the L1 norm instead of the L2 one (i.e., minimizing a sum of
# absolute values instead of the sum of squares). In turn, the uniqueness
# property of the solution is lost, but we do not need a QP-capable optimizer
# at all:

linear_solution =
    linear_parsimonious_flux_balance_analysis(model; optimizer = HiGHS.Optimizer)

#

linear_solution.fluxes

@test isapprox(linear_solution.parsimonious_objective, 518.422; atol = TEST_TOLERANCE) #src
@test isapprox( #src
    sum(abs.(values(linear_solution.fluxes))), #src
    linear_solution.parsimonious_objective, #src
    atol = TEST_TOLERANCE, #src
) #src

# ## Implementing CycleFreeFlux using parsimonious analysis
#
# CycleFreeFlux essentially defines a L1-parsimonious model which can be used
# to run a cycle-free FBA and FVA. In COBREXA, this is best done with
# [`linear_parsimonious_flux_balance_analysis`](@ref).
#
# First, let's create a constraint tree with the model, and ask for explicitly
# materializing constraints for the exchanges.

cs = flux_balance_constraints(model, interface = :identifier_prefixes)

# We will also need some existing solution of the model -- CycleFreeFlux
# algorithm uses this one as a reference for fixing the exchange reaction flux.

some_flux =
    optimized_values(cs, objective = cs.objective.value, optimizer = HiGHS.Optimizer)

# (Ideally, we should use a solving method that gives a more unique flux, but for this example a simple FBA optimum will do.)
#
# With this in hand, we can start the CycleFreeFlux workflow by placing
# constraints on exchange reactions in a linear-parsimonious model:

import ConstraintTrees as C

cs = linear_parsimonious_flux_balance_constraints(model)

cs *=
    :fixed_exchanges^C.ConstraintTree(
        k => C.Constraint(cs.fluxes[k].value, relative_tolerance_bound(0.999)(v)) for
        (k, v) in some_flux.interface.exchanges
    )

# (We purposefully made the constraints a little less strict by using
# [`relative_tolerance_bound`](@ref) -- the toy E. coli model would otherwise
# display no variability at all.)
#
# Now we can get a L1-parsimonious (thus cycle-free) solution of the model with
# the predefined exchanges:

cycle_free_flux = parsimonious_optimized_values(
    cs,
    objective = cs.objective.value,
    objective_value = some_flux.objective,
    parsimonious_objective = cs.parsimonious_objective.value,
    optimizer = HiGHS.Optimizer,
)

cycle_free_flux.fluxes

# ### CycleFreeFVA
#
# With this in hand, we can also run the cycle-free flux variability analysis
# (again with an added bit of tolerances in both the objective and parsimonious
# bounds):

cs.objective.bound = C.Between(cycle_free_flux.objective * 0.999, Inf)
cs.parsimonious_objective.bound =
    C.Between(0, cycle_free_flux.parsimonious_objective * 1.001)

var = constraints_variability(cs, cs.fluxes, optimizer = HiGHS.Optimizer)

# ### CycleFree sampling
#
# Naturally, we can also run flux sampling from the above model. To implement
# this, we follow the implementation of [`flux_sample`](@ref) --- first we
# generate the warmup:

import JuMP
warmup = vcat(
    (
        transpose(v) for (_, vs) in constraints_variability(
            cs,
            cs.fluxes,
            optimizer = HiGHS.Optimizer,
            output = (_, om) -> JuMP.value.(om[:x]),
            output_type = Vector{Float64},
        ) for v in vs
    )...,
)

# Next, we can run the sampling:

sample = sample_constraints(
    sample_chain_achr,
    cs,
    start_variables = warmup,
    seed = UInt(1234),
    output = cs.fluxes,
    n_chains = 10,
    collect_iterations = collect(100:105),
)

# The results can be observed (and usually plotted) from the sample vectors,
# such as the one for oxygen exchange:

sample.EX_o2_e
