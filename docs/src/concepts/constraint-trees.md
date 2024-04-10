
# Constraint trees

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
   (such as SBML, or lab measurements)
2. COBREXA functions are used to convert these to a constraint tree that
   properly describes the problem at hand
   - possibly, multiple types and groups of raw data can be soaked into the
     constraint tree
3. Analysis functionality of COBREXA is used to solve the system described by
   the constraitn tree, and extract useful information from the solutions.

COBREXA mainly provides functionality to make this workflow easy to use:

- [Front-end functions](#TODO) help to run the 3 above steps easily without any
  intermediate steps.
- Front-end functions call various [Front-end constraint tree builders](#TODO)
  which translate various kinds of raw data to the constraint trees
- Additional [constraint builders](#TODO) are provided to decorate the "raw"
  model representations with various additional semantics and limits
- Some functions are provided to [join the constraint trees](#TODO) via
  interfaces, simplifying e.g. the creation of community or multi-organ models,
- Finally, the [analysis function](#TODO) simulate the model in the constraint
  tree mechanistically and extract analysis results.

Technical description of the constraint tree functionality together with
examples of basic functionality and many useful utility functions is available
in [dedicated documentation of
ConstraintTrees.jl](https://cobrexa.github.io/ConstraintTrees.jl/).
