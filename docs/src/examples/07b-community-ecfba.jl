
# Copyright (c) 2021-2025, University of Luxembourg                         #src
# Copyright (c) 2021-2025, Heinrich-Heine University Duesseldorf            #src
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

# # Enzyme-constrained communities
#
# This example demonstrates a simple way to create enzyme-constrained community
# models. We use the genome-scale iML1515 E. coli model to simulate a community
# of 2 auxotrophic organisms (both of which exhibit no growth in isolation due
# to the lack of the amino acid), and observe how they can save each other by
# supplying themselves amino-acids. Such analysis easily extends to more
# auxotrophes and even communities of several different species.
#
# The simulations are, very roughly, replicating the logic of the experimental
# work by *Mee, Michael T., et al.  "Syntrophic exchange in synthetic microbial
# communities." Proceedings of the National Academy of Sciences 111.20 (2014)*.

# As usual, we start by loading packages and downloading models:

using COBREXA
import AbstractFBCModels as A
import JSONFBCModels
import ConstraintTrees as C
import HiGHS

download_model(
    "http://bigg.ucsd.edu/static/models/iML1515.json",
    "iML1515.json",
    "b0f9199f048779bb08a14dfa6c09ec56d35b8750d2f99681980d0f098355fbf5",
)

# ## Collecting data and parameters
#
# Enzyme-constrained models require parameters for protein molar masses and
# reaction turnover numbers (kcats). COBREXA supplies prepared example data for
# the iML1515 model; in this section we summarize the loading of the data into
# Julia structures from the used format. Other formats will work just as well.
#
# The loading is hidden by default for brevity:
#
#md # ```@raw html
#md # <details><summary><strong>Loading the model parameters</strong></summary>
#md # ```

import CSV

data_dir = joinpath(dirname(pathof(COBREXA)), "..", "docs", "src", "examples", "data");

e_coli_gp_mass = Dict{String,Float64}(
    x.gene_product => x.mass for
    x in CSV.File(joinpath(data_dir, "e_coli_gp_mass.tsv"), delim = '\t')
);

kcat_scale = 3600 / 1e3;
e_coli_rxn_kcat_isozyme = Dict{String,Isozyme}(
    x.reaction => Isozyme(
        kcat_forward = x.kcat * kcat_scale,
        kcat_reverse = x.kcat * kcat_scale,
        gene_product_stoichiometry = Dict(),
    ) for x in CSV.File(joinpath(data_dir, "e_coli_reaction_kcat.tsv"), delim = '\t')
);

e_coli_rxn_isozymes = Dict{String,Dict{String,Isozyme}}();
for x in CSV.File(joinpath(data_dir, "e_coli_isozyme_gp.tsv"), delim = '\t')
    haskey(e_coli_rxn_kcat_isozyme, x.reaction) || continue
    rxn = get!(e_coli_rxn_isozymes, x.reaction, Dict{String,Isozyme}())
    iso = get!(rxn, x.isozyme, deepcopy(e_coli_rxn_kcat_isozyme[x.reaction]))
    iso.gene_product_stoichiometry[x.gene_product] = x.stoichiometry
end;

#md # ```@raw html
#md # </details>
#md # ```
#
# In the end, we have gene product weight data (just like in the
# [enzyme-constrained model example](05b-enzyme-constrained-models.md)):

e_coli_gp_mass

# ... as well as isozyme data with kcats:

e_coli_rxn_isozymes

# ## Model assembly
#
# For simplicity, we will work with the "canonical" Julia-structured view of
# the iML1515:

wt_model = load_model("iML1515.json", A.CanonicalModel.Model)

# As the usual quirk, we loosen the lower bound on glucose intake that is
# required for plain FBA:
wt_model.reactions["EX_glc__D_e"].lower_bound = -1000.0;

# Additionally we allow the models isoleucine and methionine uptake:
wt_model.reactions["EX_ile__L_e"].lower_bound = -1000.0;
wt_model.reactions["EX_met__L_e"].lower_bound = -1000.0;

# ...and for good manners, we also remove the biomass annotation from the
# biomass reaction that we are not interested in:
wt_model.reactions["BIOMASS_Ec_iML1515_WT_75p37M"].annotations["sbo"] = [];

# Let's create these two knockouts-- one incapable of producing isoleucine:

ile_model = deepcopy(wt_model)
delete!(ile_model.reactions, "THRD_L");

# ...and another one without the reaction that is required for producing
# methionine:

met_model = deepcopy(wt_model)
delete!(met_model.reactions, "HSST");

# For brevity, let's make a shortcut that creates enzyme-constrained FBA system
# from the model together with a proper interface for community building:

ecfba_constraints(m, capacity) = enzyme_constrained_flux_balance_constraints(
    m,
    reaction_isozymes = e_coli_rxn_isozymes,
    gene_product_molar_masses = e_coli_gp_mass,
    interface = :sbo;
    capacity,
)

# We can now create the community by creating each model's constraint tree with
# an interface with [`enzyme_constrained_flux_balance_constraints`](@ref), and
# connecting them via [`community_flux_balance_constraints`](@ref). We have to
# pick the model abundances for cFBA, so we pick 1:1 abundance ratio. We also
# have to pick the capacities for the enzyme-constrained models (these will be
# properly diluted by the community FBA formulation), and specify that the
# community is not allowed to exchange either of our two selected amino acids
# externally (the individual models might cheat the auxotrophe community
# setting by consuming these).

community_constraints = community_flux_balance_constraints(
    [
        "ile_ko" => (ecfba_constraints(ile_model, 100.0), 0.5),
        "met_ko" => (ecfba_constraints(met_model, 100.0), 0.5),
    ],
    ["EX_ile__L_e" => 0.0, "EX_met__L_e" => 0.0],
)

# ## Simulating the community
#
# Since the community constraints created above form a completely normal
# optimization problem, we can optimize them as usual via
# [`optimized_values`](@ref); picking the `community_biomass` value as an
# objective:

res = optimized_values(
    community_constraints,
    objective = community_constraints.community_biomass.value,
    optimizer = HiGHS.Optimizer,
)

# We can observe that the community indeed grows, although not as quickly as
# the WT model normally would:

res.community_biomass
@test isapprox(res.community_biomass, 0.2304879809596633, atol = TEST_TOLERANCE) #src

# One may also observe the "global" community exchanges:
sort(collect(res.community_exchanges), by = last)

# Appropriately, we can check that the individual community members exchange
# the expected amino acids (the individual values are scaled to the individual
# members' biomasses; in the community view these values would be halved by the
# 0.5 abundances):

[
    res.met_ko.fluxes.EX_ile__L_e res.met_ko.fluxes.EX_met__L_e
    res.ile_ko.fluxes.EX_ile__L_e res.ile_ko.fluxes.EX_met__L_e
]

# We can see that isoleucine is indeed moving into the isoleucine knockout, and
# methionine into the methionine knockout. (The signs follow the usual exchange
# convention where negative values mean uptake and positive values mean
# secretion.)
#
# Finally, one might be interested in finding the optimal community composition
# for the auxotrophes. [The example on community
# building](04-community-models.md) describes a common way to find the optimal
# abundance ratios via screening.
