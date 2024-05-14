
# Logical structure of COBREXA

COBREXA uses
[ConstraintTrees.jl](https://github.com/COBREXA/ConstraintTrees.jl) for
internal representation of all metabolic modeling problems. In short,
constraint trees are "tidy" representations of the constraint-based modeling
problems, which store information about the variables' participation in named
constraints. They provide several main benefits:

- User is freed from having to allocate variables -- all variables are present
  only implicitly, are described by their "semantics" by participation in
  constraints, and are allocated automatically whenever the user connects the
  constraint trees.
- There is no need for complicated index manipulation (as with linear-algebraic
  "matrixy" model representations) nor for various identifier mangling schemes
  -- constraint trees provide named interface for everything, and identifiers
  can be organized into directories to prevent name clashes in various
  multi-compartment and community models.
- Contrary to the fixed model representations (such as SBML or JSON models),
  ConstraintTrees do not possess a semantic of a "single flux-based model" and
  are thus infinitely extensible, allowing easy creation, manipulation and
  storage of even very complicated constraint systems.

With ConstraintTrees, the typical workflow in COBREXA is as follows:

1. "Raw" data and base model data are loaded from semantically organized models
   (such as SBML, or lab measurements in CSV or other tabular format)
2. COBREXA functions are used to convert these to a constraint tree that
   properly describes the problem at hand
   - possibly, multiple types and groups of raw data can be soaked into the
     constraint tree
3. Analysis functionality of COBREXA is used to solve the system described by
   the constraitn tree, and extract useful information from the solutions.

COBREXA mainly provides functionality to make this workflow easy to use for
many various purposes:

- **Front-end functions** help to run the 3 above steps easily without any
  intermediate steps. These include:
  - [`community_flux_balance_analysis`](@ref)
  - [`enzyme_constrained_flux_balance_analysis`](@ref)
  - [`flux_balance_analysis`](@ref)
  - [`flux_sample`](@ref)
  - [`gene_knockouts`](@ref)
  - [`linear_metabolic_adjustment_minimization_analysis`](@ref)
  - [`loopless_flux_balance_analysis`](@ref)
  - [`max_min_driving_force_analysis`](@ref)
  - [`metabolic_adjustment_minimization_analysis`](@ref)
  - [`objective_production_envelope`](@ref)
  - [`simplified_enzyme_constrained_flux_balance_analysis`](@ref)
- Front-end functions call various **Front-end constraint tree builders** which
  translate various kinds of raw data to the constraint trees, such as:
  - [`community_flux_balance_constraints`](@ref)
  - [`enzyme_constrained_flux_balance_constraints`](@ref)
  - [`flux_balance_constraints`](@ref)
  - [`flux_variability_analysis`](@ref)
  - [`gene_knockout_constraints`](@ref)
  - [`linear_metabolic_adjustment_minimization_constraints`](@ref)
  - [`log_concentration_constraints`](@ref)
  - [`loopless_flux_balance_constraints`](@ref)
  - [`max_min_driving_force_constraints`](@ref)
  - [`metabolic_adjustment_minimization_constraints`](@ref)
  - [`simplified_enzyme_constrained_flux_balance_constraints`](@ref)
- Additional **constraint builders** are provided to decorate the "raw" model
  representations with various additional semantics and limits:
  - [`all_equal_constraints`](@ref), [`greater_or_equal_constraint`](@ref) and
    similar ones
  - [`knockout_constraints`](@ref)
  - [`loopless_constraints`](@ref)
  - [`scale_bounds`](@ref)
  - [`scale_constraints`](@ref)
  - [`sign_split_variables`](@ref) and [`sign_split_constraints`](@ref)
    [`squared_sum_error_value`](@ref)
  - [`sum_value`](@ref), [`squared_sum_value`](@ref) and
- Some functions are provided to **join the constraint trees** via interfaces,
  simplifying e.g. the creation of community or multi-organ models,
  - [`flux_balance_constraints`](@ref) can automatically generate interfaces
    suitable for community-style and multi-compartment-style metabolic
    modeling from the annotations in the FBC models
  - [`interface_constraints`](@ref) joins the "modules" with prepared
    interfaces together
- Finally, the **analysis functions** simulate the model in the constraint
  tree mechanistically and extract analysis results:
  - [`constraints_objective_envelope`](@ref)
  - [`parsimonious_optimized_values`](@ref)
  - [`sample_constraints`](@ref)
  - [`sample_constraint_variables`](@ref)
  - [`screen_optimization_model`](@ref)
  - [`screen`](@ref)
  - [`optimized_values`](@ref)
  - [`constraints_variability`](@ref)

!!! tip "Exploring and customizing the frontend analysis functions"
    To know which builder function is used to create or modify some kind of
    constraint tree in COBREXA, use the "link to source code" feature in the
    frontend function's individual documentation. The source code of front-end
    functions is written to be as easily re-usable as possible -- one can
    simply copy-paste it into the program, and immediately start building
    specialized and customized front-end functions.

Technical description of the constraint tree functionality, together with
examples of basic functionality and many useful utility functions is available
in [dedicated documentation of
ConstraintTrees.jl](https://cobrexa.github.io/ConstraintTrees.jl/).
