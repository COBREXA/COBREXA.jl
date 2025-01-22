
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

# # Simplified resource balance analysis
#
# Resource balance analysis (RBA) models are extensions of enzyme constrained
# models that additionally incorporate other cellular mechanisms, such as
# translation, transcription and replication. This requires much more
# mechanistic knowledge about the processes, but also dramatically improves the
# predictive capability of the model.
#
# Here we demonstrate the approach for building such extensions with COBREXA,
# over a demonstrational simplified RBA model that accounts for the major
# translation costs (synthesis of proteins and ribosomes).
#
# For comprehensiveness, we use the full genome-scale model of E. coli
# (iML1515):

using COBREXA

download_model(
    "http://bigg.ucsd.edu/static/models/iML1515.json",
    "iML1515.json",
    "b0f9199f048779bb08a14dfa6c09ec56d35b8750d2f99681980d0f098355fbf5",
)

# We use several packages as usual; additionally we import `ConstraintTrees`
# for later modifications.

import AbstractFBCModels as A
import JSONFBCModels
import HiGHS
import ConstraintTrees as C

# ## Data and parameters for the RBA model
#
# For the purposes of this example, COBREXA comes with example data for the
# whole iML1515 model, aggregated from several publications and databases. In
# this section we simply load the data into suitable Julia structures. Other
# data formats may work just as well.
#
# The loading is hidden by default for brevity:
#
#md # ```@raw html
#md # <details><summary><strong>Loading the RBA model parameters</strong></summary>
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

e_coli_gp_aas = Dict{String,Dict{Symbol,Int}}(
    begin
        d = Dict(keys(x) .=> values(x))
        gp = d[:gene_product]
        delete!(d, :gene_product)
        gp => d
    end for x in CSV.File(joinpath(data_dir, "e_coli_gp_aa.tsv"), delim = '\t')
);

amino_acids = Set(Symbol(aa) for (k, v) in e_coli_gp_aas for (aa, _) in v);

#md # ```@raw html
#md # </details>
#md # ```
#
# In the end, we have gene product weight data (just like in the [enzyme-constrained model example](05b-enzyme-constrained-models.md)):

e_coli_gp_mass

# ... as well as isozyme data with kcats:

e_coli_rxn_isozymes

# We additionally need a list of amino acids in the model:

amino_acids

# ...together with a list of how much amino acids there are in which gene product:

e_coli_gp_aas

# To make the RBA problem working, we also need to assume some constant parameters (many of such can be found via https://bionumbers.hms.harvard.edu):

atp_polymerization_cost = 0.042;
protein_polymerization_atp_per_gDW = 12.0;
ribosome_speed_aa_per_hour = 12.0 * 3600;
ribosome_molar_mass = 2700.0;

# We also need a stoichiometry for "energy consumption" reaction, which we will
# use to simulate the energy cost of translation:
energy_stoichiometry = Dict(:atp_c => -1, :h2o_c => -1, :adp_c => 1, :h_c => 1, :pi_c => 1)

# ## Model assembly
# ### Enzyme-constrained base model
#
# First, we load the model in a format that is suitable for doing small
# changes:
model = load_model("iML1515.json")

# We will require some access to the stoichiometry of the biomass reaction (in
# essence, we copy a part of it, but replace the part that uses amino acids as
# a building material, and slightly enhance the energy consumption part). So we
# save it here:
biomass = Dict(
    Symbol(k) => v for
    (k, v) in A.reaction_stoichiometry(model, "BIOMASS_Ec_iML1515_core_75p37M")
)

# We can create the enzyme-constrained model for iML1515. This will be extended
# later.
ec_constraints = enzyme_constrained_flux_balance_constraints(
    model;
    reaction_isozymes = e_coli_rxn_isozymes,
    gene_product_molar_masses = e_coli_gp_mass,
    capacity = 550.0,
)

# Before we continue, we apply a small quirk to remove an artificial limit on
# the glucose intake. (The limit is required to prevent "infinite" growth in
# simplistic FBA-style analysis; in our case the enzyme capacity serves as a
# sufficient and more realistic limiter. We have to unblock both the
# bidirectional reaction and the "reversed" view, since both carry the bound.)
ec_constraints.fluxes.EX_glc__D_e.bound.lower = -1000
ec_constraints.fluxes_reverse.EX_glc__D_e.bound.upper = 1000

# To avoid the model from growing in unexpected modes, we will constraint the
# original biomass reaction to zero:
ec_constraints.fluxes.BIOMASS_Ec_iML1515_core_75p37M.bound = C.EqualTo(0)
ec_constraints.fluxes.BIOMASS_Ec_iML1515_WT_75p37M.bound = C.EqualTo(0)

# ### RBA translation machinery

# A common issue with RBA formulations is that the biomass-based growth formula
# depends on a determined optimal composition of the enzyme pool and actual
# production of metabolites; which makes the underlying constrained problem
# quadratic.
#
# A common way to dodge the need for quadratic solvers is to solve the problem
# for a fixed growth rate, which we is the approach that we choose here.
# Alternatively, one might state the full quadratic problem and solve it, with
# some performance cost stemming from use of QP solvers.
#
# Let's first make a utility function that prepares the connection to the
# metabolite pool, and adds several useful variables atop a given
# enzyme-constrained model:
function with_translation_variables(ec_constraints::C.ConstraintTree)
    #+
    # Create a "resource pool" and connect it to the stoichiometry of the
    # intracellular (and other) metabolites. (This effectively creates new
    # exchange reactions in `ec_constraints`.)
    resources = C.variables(
        keys = Symbol.(
            collect(union(amino_acids, keys(energy_stoichiometry), keys(biomass)))
        ),
    )
    (cs, rs) =
        inject_interface(ec_constraints, :flux_stoichiometry^resources, multiplier = -1)
    #+
    # Also add a single new variable for the production of ribosomes by
    # ribosomes.
    return cs * :resources^rs.flux_stoichiometry +
           :ribosome_production^C.variable(; bound = (0, Inf))
end

# To make the construction nicer, we'll make a helper for summing up
# constraint-tree values:
sum_values(x...) = C.sum(x..., init = zero(C.LinearValue))

# ...and another helper for adding values in constraint trees together:
add_trees(ts...) =
    C.preduce(ts, init = C.ConstraintTree()) do t1, t2
        z(::Missing) = zero(C.LinearValue)
        z(x) = C.value(x)
        return C.merge(t1, t2) do c1, c2
            C.Constraint(z(c1) + z(c2))
        end
    end

# Since we have to solve the problem for multiple growth rates to be able to
# scan for optimum, we will wrap the growth-dependent part in a reusable
# function:
function translation_constraints(
    resources::C.ConstraintTree,
    ribos_required_for_ribos::C.Constraint,
    gene_product_amounts::C.ConstraintTree,
    growth::Float64,
)
    #+
    # First we can calculate how much amino acids we need to build the expected
    # amount of gene products:
    aas_required_for_gps = C.imap(gene_product_amounts) do (i,), gp
        C.ConstraintTree(
            aa => C.Constraint(gp.value * v * growth * 0.001) for
            (aa, v) in e_coli_gp_aas[String(i)]
        )
    end
    #+
    # This allows us to calculate how much ribosome we have to produce to make
    # the production of all of the above enzymes possible:
    ribo_required_for_gps = C.ConstraintTree(
        i => C.Constraint(
            sum_values(aa.value for (_, aa) in aas) / ribosome_speed_aa_per_hour,
        ) for (i, aas) in aas_required_for_gps
    )
    #+
    # Now we know the total amount of ribosome to produce (both for the above
    # protein production and for production of ribosomes itself) so we can see
    # how much amino acids in total are required for production of ribosomes:
    total_ribos_required =
        ribos_required_for_ribos.value +
        sum_values(c.value for (_, c) in ribo_required_for_gps)
    aas_required_for_ribos = C.ConstraintTree(
        aa => C.Constraint(total_ribos_required * v) for
        (aa, v) in e_coli_gp_aas["ribosome"]
    )
    #+
    # Now we solve a "rocket equation" -- the ribosomes need to produce both
    # the protein-producing ribosomes and themselves, so we add a constraint
    # that ensures there's enough ribosomes for both.
    ribosome_balance_constraint = equal_value_constraint(
        sum_values(aa.value for (_, aa) in aas_required_for_ribos),
        ribos_required_for_ribos.value * ribosome_speed_aa_per_hour,
    )
    #+
    # With the AA requirements solved, we can estimate how much energy we need
    # for polymerization of the proteins and ribosomes:
    energy_required =
        atp_polymerization_cost *
        sum_values(aa.value for (_, aas) in aas_required_for_gps for (_, aa) in aas) +
        sum_values(v.value for (_, v) in aas_required_for_ribos)
    #+
    # With all that in hand, we can put together the final resource consumption:
    resource_consumption = add_trees(
        C.values(aas_required_for_gps)...,
        aas_required_for_ribos,
        C.map(stoi -> -stoi * energy_required, C.Tree{Int}(energy_stoichiometry), C.Value),
    )
    #+
    # ...and make a stoichiometry out of that, with exact cases for amino acids
    # (these are completely replaced in the original biomass), energy
    # metabolites (these are partially re-used from the original biomass, but
    # with an adjustment that tries to remove the polymerization cost portion
    # in the original model), and everything other scaled for growth:
    resource_stoichiometry = C.imap(resources) do (resource,), input
        if resource in amino_acids # AA case
            equal_value_constraint(resource_consumption[resource], input)
        elseif resource in keys(energy_stoichiometry) # energy case
            equal_value_constraint(
                -growth * (
                    biomass[resource] -
                    energy_stoichiometry[resource] * protein_polymerization_atp_per_gDW
                ) -
                resource_consumption[resource].value * energy_stoichiometry[resource],
                input,
            )
        else # everything else
            C.Constraint(input.value, -growth * biomass[resource])
        end
    end
    #+
    # Finally, let's wrap all the constraints and some useful derived helper
    # values in one big tree:
    return C.ConstraintTree(
        :resource_stoichiometry => resource_stoichiometry,
        :ribosome_balance => ribosome_balance_constraint,
        :gene_product_production => ribo_required_for_gps,
        :total_ribosome_mass => C.Constraint(
            ribosome_molar_mass * (
                sum_values(v.value for (_, v) in ribo_required_for_gps) +
                ribos_required_for_ribos.value
            ),
        ),
        :amino_acid_use => (aas_required_for_gps * :ribosome^aas_required_for_ribos),
        :polymerization_energy => C.Constraint(energy_required),
        :translation_resource_consumption => resource_consumption,
    )
end

# ## Running the resource-balanced simulation

# With the above functions, assembling a resource-balanced model amounts to
# adding new variables and connecting them with the rest of the
# enzyme-constrained model. We assemble a model for growth value of 0.6
# gDW/gDWh:
rb_constraints = with_translation_variables(ec_constraints)
rb_constraints *= translation_constraints(
    rb_constraints.resources,
    rb_constraints.ribosome_production,
    rb_constraints.gene_product_amounts,
    0.9,
)

# The model may be slightly under-constrained for less-than-extreme values of
# growth; to obtain a realistic solution for we can ask the solver to minimize
# the mass of used resources. Accordingly, we re-constraint the total mass of
# the model:
rb_constraints.gene_product_capacity.total_capacity.bound = nothing
rb_constraints.total_capacity = C.Constraint(
    rb_constraints.gene_product_capacity.total_capacity.value +
    rb_constraints.total_ribosome_mass.value,
    (0.0, 550.0),
)

# We can optimize the model now, minimizing the mass:
res = optimized_values(
    rb_constraints,
    objective = rb_constraints.total_capacity.value,
    sense = Minimal,
    optimizer = HiGHS.Optimizer,
)

# The model can be used to observe various interesting effects. For example,
# how much building material is required to reach such growth?
res.total_capacity
@test isapprox(res.total_capacity, 538.778495150974, atol = TEST_TOLERANCE) #src

# What is the resource composition used for building up biomass?
sort(collect(res.resources), by = last)

# How much of that comes into (and out of) translation?
res.translation_resource_consumption

# How much arginine is used to build the ribosomes?
res.amino_acid_use.ribosome.arg__L_c

# This solution is not necessarily optimal though. To find an optimal growth,
# one may use e.g. [`screen`](@ref) to run the same simulation over many growth
# values, and pick the largest feasible growth.
