
# Optimizer compatibility and settings

COBREXA uses [JuMP.jl](https://github.com/jump-dev/JuMP.jl) as an abstraction
layer above the Julia optimization ecosystem and
[MathOptInterface.jl](https://github.com/jump-dev/MathOptInterface.jl).
Problems represented in COBREXA are typically not implemented in a
optimizer-specific manner, and it is very easy to switch between different optimizers
and various configurations to try which ones are able to represent the model
best. Generally, you can install any optimizer from [JuMP's
list](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) and work
with it.

For example, to install SCIP, we request it from Julia packages:
```
] add SCIP
```
...after which we can use `SCIP.Optimizer` object as a parameter to anything
that requires the optimizer:
```julia
using SCIP
flux_balance_analysis(some_model, optimizer = SCIP.Optimizer)
```

## Optimizer choice considerations

The desired optimizer is typically a subject to various constraints. Most
notably, the optimizer must be capable of processing the optimization problem
at hand. The choices reduce as follows:

- For analyses that require quadratic objectives (typically the
  [L2-parsimonious analyses](examples/03b-parsimonious-flux-balance.md)), an
  optimizer with QP support is required.
- [Gap-filling](@ref gap_filling_analysis), [medium optimization](@ref
  medium_optimization_analysis) and [loopless FBA](@ref
  loopless_flux_balance_analysis) typically require MILP support in the
  optimizer.
- In case you require numerically stable solutions, an IPM optimizer may be
  viable.

Notably, some of the "larger" optimizers may support all of these options
(these include SCIP, HiGHS, and the commercial Gurobi and CPLEX) at the slight
cost of complexity and potential slowness in general cases.

## Optimizer parameters

Most optimizer may be tuned by attributes, set via `settings` argument in most
optimizing functions, and via optimizer attributes functionality of JuMP. In
COBREXA, the attributes can be modified via [`set_optimizer_attribute`](@ref).

Unfortunately, the settings of individual optimizers differ quite
significantly, so there is no good universal default and one has to consult the
documentation of the given optimizer to find the actual names and expected
values of attributes. On the other hand, following some simple guidelines may
push most optimizers into behaving nearly optimally. We list the common
considerations here:

- Tolerance settings do not have to be too strict in most cases; tolerance of
  `1e-5` is typically good for working with genome-scale models.
  - In some cases, the tolerance has to be made stricter; e.g., gap-filling may
    "enable" a little amount of flux through an universal reaction that slips
    through the tolerance even without setting the corresponding binary
    variable to `1`, which in turn yields invalid solution. **Tolerance
    settings in potentially numerically unstable MILP programs should always be
    validated**.
  - In large models (e.g., the [enzyme-constrained
    communities](examples/07b-community-ecfba.md)), individual
    in-tolerance errors may accumulate and harshly impact the precision of the
    solution. Apart from tightening the tolerance, some models may also require
    rescaling in order to be more numerically stable.
- Some analysis types (notably [FVA](@ref flux_variability_analysis), [warm-up
  computation for samplers](@ref flux_sample), but also [gene knockouts](@ref
  gene_knockouts) and [production envelopes](@ref
  objective_production_envelope)) benefit from the ability of the solver to
  quickly restart computation with only small modifications from the previous
  known solution (usually this is called a "hot start"). Especially for FVA,
  this brings great performance benefits. For such cases, one often has to
  force the solver to use a seemingly sub-optimal solving strategy that is
  sufficiently simple to allow the fast restarts.
  - As a general rule of thumb, any kind of dual-solving method should be
    disabled for FVA and sampler warm-up, because the change of the
    optimization objective between algorithm iterations implies recomputation
    of the whole dual program. In particular **for FVA, we obtained best results
    with the "primal simplex" solving methods**.
  - For computations such as gene knockout analysis, the primal simplex is
    normally changed in each algorithm iteration, forcing a recomputation
    again. We obtained interesting results by using IPM solvers in such cases
    -- the IPM solution often stays valid between different knockouts, speeding
    up the whole computation substantially.
- When solving linear programs in parallel (e.g., on a cluster or HPC), the
  additional parallelism introduced in the solver itself may cause scheduling
  issues, resource starvation and eventually general slowness. In such cases,
  reducing the number of solver threads to a single one may actually help the
  performance.

For example, a good FVA-capable settings for HiGHS (gathered from [its
documentation](https://ergo-code.github.io/HiGHS/dev/options/definitions/)) may
look like this:

```julia
import HiGHS
optimizer = HiGHS.Optimizer
settings = [
    set_optimizer_attribute("solver", "simplex"),
    set_optimizer_attribute("simplex_strategy", 4), # stands for "primal simplex"
    set_optimizer_attribute("parallel", "off"),
]

# ...

flux_variability_analysis(my_model; optimizer, settings)
```
