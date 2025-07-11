
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

# # Enzyme constrained models
# Enzyme constrained metabolic models include the effect of enzyme kinetics (v =
# k * e) and a protein capacity limitation (∑e = Etotal) on conventional mass
# balance (FBA) models.

using COBREXA

# Here we will construct an enzyme constrained variant of the *E. coli* "core"
# model. We will need the model, which we can download if it is not already present.

download_model(
    "http://bigg.ucsd.edu/static/models/e_coli_core.json",
    "e_coli_core.json",
    "7bedec10576cfe935b19218dc881f3fb14f890a1871448fc19a9b4ee15b448d8",
)

# Additionally to COBREXA and the model format package, we will need a solver
# -- let's use HiGHS here:

import AbstractFBCModels as A
import ConstraintTrees as C
import JSONFBCModels
import HiGHS

model = load_model("e_coli_core.json")

# Enzyme constrained models require parameters that are usually not used by
# conventional constraint based models. These include reaction specific turnover
# numbers, molar masses of enzymes, and protein capacity bounds.

# ## Reaction turnover numbers

# Enzyme constrained models require reaction turnover numbers, which are often
# isozyme-specific. Many machine learning tools, or experimental data sets, can
# be used to estimate these parameters.

#md # ```@raw html
#md # <details><summary><strong>Data for reaction turnover numbers</strong></summary>
#md # ```
# This data is taken from: *Heckmann, David, et al. "Machine learning applied
# to enzyme turnover numbers reveals protein structural correlates and improves
# metabolic models." Nature communications 9.1 (2018): 1-10.*
const ecoli_core_reaction_kcats = Dict(
    "ACALD" => 568.11,
    "PTAr" => 1171.97,
    "ALCD2x" => 75.95,
    "PDH" => 529.76,
    "MALt2_2" => 234.03,
    "CS" => 113.29,
    "PGM" => 681.4,
    "TKT1" => 311.16,
    "ACONTa" => 191.02,
    "GLNS" => 89.83,
    "ICL" => 17.45,
    "FBA" => 373.42,
    "FORt2" => 233.93,
    "G6PDH2r" => 589.37,
    "AKGDH" => 264.48,
    "TKT2" => 467.42,
    "FRD7" => 90.20,
    "SUCOAS" => 18.49,
    "ICDHyr" => 39.62,
    "AKGt2r" => 234.99,
    "GLUSy" => 33.26,
    "TPI" => 698.30,
    "FORt" => 234.38,
    "ACONTb" => 159.74,
    "GLNabc" => 233.80,
    "RPE" => 1772.485,
    "ACKr" => 554.61,
    "THD2" => 24.73,
    "PFL" => 96.56,
    "RPI" => 51.77,
    "D_LACt2" => 233.51,
    "TALA" => 109.05,
    "PPCK" => 218.42,
    "PGL" => 2120.42,
    "NADTRHD" => 186.99,
    "PGK" => 57.64,
    "LDH_D" => 31.11,
    "ME1" => 487.01,
    "PIt2r" => 233.86,
    "ATPS4r" => 71.42,
    "GLCpts" => 233.90,
    "GLUDy" => 105.32,
    "CYTBD" => 153.18,
    "FUMt2_2" => 234.37,
    "FRUpts2" => 234.19,
    "GAPD" => 128.76,
    "PPC" => 165.52,
    "NADH16" => 971.74,
    "PFK" => 1000.46,
    "MDH" => 25.93,
    "PGI" => 468.11,
    "ME2" => 443.09,
    "GND" => 240.12,
    "SUCCt2_2" => 234.18,
    "GLUN" => 44.76,
    "ADK1" => 111.64,
    "SUCDi" => 680.31,
    "ENO" => 209.35,
    "MALS" => 252.75,
    "GLUt2r" => 234.22,
    "PPS" => 706.14,
    "FUM" => 1576.83,
)
#md # ```@raw html
#md # </details>
#md # ```

# We have these here:

ecoli_core_reaction_kcats # units = 1/s

# Each reaction in a constraint-based model usually has gene reaction rules
# associated with it. These typically take the form of, possibly multiple,
# isozymes that can catalyze a reaction. A turnover number needs to be assigned
# to each isozyme, as shown below. Additionally, some enzymes are composed of
# multiple subunits, which differ in subunit stoichiometry. This also needs to
# be accounted for. Assuming a stoichiometry of 1 for everything tends to work
# just right OK if there is no better information available.

reaction_isozymes = Dict{String,Dict{String,Isozyme}}() # a mapping from reaction IDs to isozyme IDs to isozyme structs.
for rid in A.reactions(model)
    grrs = A.reaction_gene_association_dnf(model, rid)
    isnothing(grrs) && continue # skip if no grr available
    haskey(ecoli_core_reaction_kcats, rid) || continue # skip if no kcat data available
    for (i, grr) in enumerate(grrs)
        d = get!(reaction_isozymes, rid, Dict{String,Isozyme}())
        d["isozyme_"*string(i)] = Isozyme( # each isozyme gets a unique name
            gene_product_stoichiometry = Dict(grr .=> fill(1.0, size(grr))), # assume subunit stoichiometry of 1 for all isozymes
            kcat_forward = ecoli_core_reaction_kcats[rid] * 3.6, # forward reaction turnover number units = 1/h
            kcat_reverse = ecoli_core_reaction_kcats[rid] * 3.6, # reverse reaction turnover number units = 1/h
        )
    end
end

#md # !!! tip "Turnover number units"
#md #     Take care with the units of the turnover numbers. In literature they are usually reported in 1/s. However, flux units are typically mmol/gDW/h, suggesting to rescale the turnover numbers to 1/h in order to use the conventional flux units.

#md # !!! tip "Missing turnover numbers"
#md #     In case the turnover numbers are missing for a single direction of the reaction catalysis, we can use `nothing` instead of the turnover number. This will prevent generation of any enzyme capacity constraints for the corresponding reaction direction.

# ## Enzyme molar masses

# We also require the mass of each enzyme, to properly weight the contribution
# of each flux/isozyme in the capacity bound(s). These data can typically be
# found in uniprot.

#md # ```@raw html
#md # <details><summary><strong>Gene product masses</strong></summary>
#md # ```
# This data is downloaded from Uniprot for E. coli K12, gene mass in kDa. To
# obtain these data manually, go to [Uniprot](https://www.uniprot.org/) and
# search using these terms: `reviewed:yes AND organism:"Escherichia coli
# (strain K12) [83333]"`.
const ecoli_core_gene_product_masses = Dict(
    "b4301" => 23.214,
    "b1602" => 48.723,
    "b4154" => 65.972,
    "b3236" => 32.337,
    "b1621" => 56.627,
    "b1779" => 35.532,
    "b3951" => 85.96,
    "b1676" => 50.729,
    "b3114" => 85.936,
    "b1241" => 96.127,
    "b2276" => 52.044,
    "b1761" => 48.581,
    "b3925" => 35.852,
    "b3493" => 53.389,
    "b3733" => 31.577,
    "b2926" => 41.118,
    "b0979" => 42.424,
    "b4015" => 47.522,
    "b2296" => 43.29,
    "b4232" => 36.834,
    "b3732" => 50.325,
    "b2282" => 36.219,
    "b2283" => 100.299,
    "b0451" => 44.515,
    "b2463" => 82.417,
    "b0734" => 42.453,
    "b3738" => 30.303,
    "b3386" => 24.554,
    "b3603" => 59.168,
    "b2416" => 63.562,
    "b0729" => 29.777,
    "b0767" => 36.308,
    "b3734" => 55.222,
    "b4122" => 60.105,
    "b2987" => 53.809,
    "b2579" => 14.284,
    "b0809" => 26.731,
    "b1524" => 33.516,
    "b3612" => 56.194,
    "b3735" => 19.332,
    "b3731" => 15.068,
    "b1817" => 35.048,
    "b1603" => 54.623,
    "b1773" => 30.81,
    "b4090" => 16.073,
    "b0114" => 99.668,
    "b3962" => 51.56,
    "b2464" => 35.659,
    "b2976" => 80.489,
    "b1818" => 27.636,
    "b2285" => 18.59,
    "b1702" => 87.435,
    "b1849" => 42.434,
    "b1812" => 50.97,
    "b0902" => 28.204,
    "b3403" => 59.643,
    "b1612" => 60.299,
    "b1854" => 51.357,
    "b0811" => 27.19,
    "b0721" => 14.299,
    "b2914" => 22.86,
    "b1297" => 53.177,
    "b0723" => 64.422,
    "b3919" => 26.972,
    "b3115" => 43.384,
    "b4077" => 47.159,
    "b3528" => 45.436,
    "b0351" => 33.442,
    "b2029" => 51.481,
    "b1819" => 30.955,
    "b0728" => 41.393,
    "b2935" => 72.212,
    "b2415" => 9.119,
    "b0727" => 44.011,
    "b0116" => 50.688,
    "b0485" => 32.903,
    "b3736" => 17.264,
    "b0008" => 35.219,
    "b3212" => 163.297,
    "b3870" => 51.904,
    "b4014" => 60.274,
    "b2280" => 19.875,
    "b2133" => 64.612,
    "b2278" => 66.438,
    "b0118" => 93.498,
    "b2288" => 16.457,
    "b3739" => 13.632,
    "b3916" => 34.842,
    "b3952" => 32.43,
    "b2925" => 39.147,
    "b2465" => 73.043,
    "b2297" => 77.172,
    "b2417" => 18.251,
    "b4395" => 24.065,
    "b3956" => 99.063,
    "b0722" => 12.868,
    "b2779" => 45.655,
    "b0115" => 66.096,
    "b0733" => 58.205,
    "b1478" => 35.38,
    "b2492" => 30.565,
    "b0724" => 26.77,
    "b0755" => 28.556,
    "b1136" => 45.757,
    "b2286" => 68.236,
    "b0978" => 57.92,
    "b1852" => 55.704,
    "b2281" => 20.538,
    "b2587" => 47.052,
    "b2458" => 36.067,
    "b0904" => 30.991,
    "b1101" => 50.677,
    "b0875" => 23.703,
    "b3213" => 52.015,
    "b2975" => 58.92,
    "b0720" => 48.015,
    "b0903" => 85.357,
    "b1723" => 32.456,
    "b2097" => 38.109,
    "b3737" => 8.256,
    "b0810" => 24.364,
    "b4025" => 61.53,
    "b1380" => 36.535,
    "b0356" => 39.359,
    "b2277" => 56.525,
    "b1276" => 97.677,
    "b4152" => 15.015,
    "b1479" => 63.197,
    "b4153" => 27.123,
    "b4151" => 13.107,
    "b2287" => 25.056,
    "b0474" => 23.586,
    "b2284" => 49.292,
    "b1611" => 50.489,
    "b0726" => 105.062,
    "b2279" => 10.845,
    "s0001" => 0.0,
)
#md # ```@raw html
#md # </details>
#md # ```

# We have the molar masses here:

ecoli_core_gene_product_masses # unit kDa = kg/mol

#md # !!! tip "Molar mass units"
#md #     Just as with the turnover numbers, take extreme care about the units of the molar masses. In literature they are usually reported in Da or kDa (g/mol). However, as noted above, flux units are typically mmol/gDW/h. Since the enzyme kinetic equation is `v = k * e` (where `k` is the turnover number) it suggests that the enzyme variable will have units of mmol/gDW. The molar masses come into play when setting the capacity limitations, e.g. usually a sum over all enzymes weighted by their molar masses as `e * M`. Thus, if the capacity limitation has units of g/gDW, then the molar masses must have units of g/mmol (i.e., kDa).

# ## Capacity limitation

# The capacity limitation usually denotes an upper bound of protein available to
# the cell. Multiple capacity bounds can be used (cytosol, membrane, etc).

total_enzyme_capacity = 50.0 # mg of enzyme/gDW

# ## Running a basic enzyme constrained model

# With all the parameters specified, we can directly use the enzyme constrained
# convenience function to run enzyme constrained FBA in one shot:

ec_solution = enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses = ecoli_core_gene_product_masses,
    capacity = total_enzyme_capacity,
    optimizer = HiGHS.Optimizer,
)

# We can notice that the objective function is a little lower than with
# unconstrained E. coli core:

ec_solution.objective

# One can also observe many interesting thing, e.g. the amount of gene product
# material required for the system to run. Importantly, the units of these
# values depend on the units used to set the turnover numbers and protein molar
# masses.

ec_solution.gene_product_amounts

# The total amount of required gene product mass is, by default, present as
# `total_capacity`:

ec_solution.gene_product_capacity

#src these values should be unique (glucose transporter is the only way to get carbon into the system)
@test isapprox(ec_solution.objective, 0.706993382849705, atol = TEST_TOLERANCE) #src
@test isapprox( #src
    ec_solution.gene_product_capacity.total_capacity, #src
    50.0, #src
    atol = TEST_TOLERANCE, #src
) #src
@test isapprox(ec_solution.fluxes.EX_glc__D_e, -10, atol = TEST_TOLERANCE) #src
@test isapprox( #src
    ec_solution.gene_product_amounts.b2417, #src
    0.011875920383431717, #src
    atol = TEST_TOLERANCE, #src
) #src

# ## Simplified models
#
# Because most active reactions typically only use a single isozyme, we may also
# use a simplified representation of the problem where this fact is reflected,
# saving the variable allocation for the isozymes.
#
# [`simplified_enzyme_constrained_flux_balance_analysis`](@ref) takes similar
# arguments as the [`enzyme_constrained_flux_balance_analysis`](@ref), but
# automatically chooses the "fastest" reaction isozyme for each reaction
# direction and builds the model with that.
#
# We additionally show how to specify more complex bounds for the capacities;
# in this case we list all fluxes for which we have isozyme data and add a
# realistic lower limit for the capacity.

minimum_enzyme_capacity = 20.0

simplified_ec_solution = simplified_enzyme_constrained_flux_balance_analysis(
    model;
    reaction_isozymes,
    gene_product_molar_masses = ecoli_core_gene_product_masses,
    capacity = Dict(
        :total => (
            Symbol.(keys(reaction_isozymes)),
            C.Between(minimum_enzyme_capacity, total_enzyme_capacity),
        ),
    ),
    optimizer = HiGHS.Optimizer,
)

# In this case, the result is the same as with the full analysis:

simplified_ec_solution.capacity_limits.total

# Gene product amounts are not present in the model but are reconstructed
# nevertheless (they are uniquely determined by the flux):

simplified_ec_solution.gene_product_amounts

@test isapprox( #src
    ec_solution.objective, #src
    simplified_ec_solution.objective, #src
    atol = TEST_TOLERANCE, #src
) #src

# ## Variability analysis with enzyme constraints
#
# Enzyme-constrained variability analysis can be executed on a model by
# combining [`enzyme_constrained_flux_balance_constraints`](@ref) (or
# [`simplified_enzyme_constrained_flux_balance_constraints`](@ref)) with
# [`constraints_variability`](@ref) (or any other analysis function):

ec_system = enzyme_constrained_flux_balance_constraints(
    model;
    reaction_isozymes,
    gene_product_molar_masses = ecoli_core_gene_product_masses,
    capacity = total_enzyme_capacity,
)

# Here, we can do the FVA "manually", first solving the system:

ec_optimum = optimized_values(
    ec_system,
    output = ec_system.objective,
    objective = ec_system.objective.value,
    optimizer = HiGHS.Optimizer,
)

# ...then creating a system constrained to near-optimal growth:

ec_system.objective.bound = C.Between(0.99 * ec_optimum, Inf)

# ...and finally, finding the extremes of the near-optimal part of the feasible
# space:

ec_variabilities =
    constraints_variability(ec_system, ec_system, optimizer = HiGHS.Optimizer)

# By default, the result computes variabilities of all possible values in the
# model. (I.e., it also computes variabilities for the variable combinations
# that are present in the tree!) As usual, the results can be observed in the
# original constraint tree structure, giving us the variabilities for reaction
# fluxes:

ec_variabilities.fluxes

# ...as well as for gene product requirements:

ec_variabilities.gene_product_amounts

# ...and for the individual directional isozymes:

ec_variabilities.isozyme_forward_amounts.PGM

# If we do not need to compute all these values, it is often more efficient to
# only ask for the part of the output that is required:

ec_gp_amount_variabilities = constraints_variability(
    ec_system,
    ec_system.gene_product_amounts,
    optimizer = HiGHS.Optimizer,
)

@test isapprox(ec_gp_amount_variabilities.b0008[1], 0, atol = TEST_TOLERANCE) #src
@test isapprox( #src
    ec_gp_amount_variabilities.b0008[2], #src
    0.020956009969910837, #src
    atol = TEST_TOLERANCE, #src
) #src
