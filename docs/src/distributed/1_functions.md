
# Parallel processing overview

Distributed processing in Julia is represented mainly by the package
[`Distributed.jl`](https://docs.julialang.org/en/v1/stdlib/Distributed/).

`COBREXA.jl` is able to utilize this existing system to almost transparently
run the large parallelizable analyses on multiple CPU cores and multiple
computers connected through the network. Ultimately, the approach scales to
thousands of computing nodes in large HPC facilities.

You may run your analyses in parallel to gain speed-ups. The usual workflow in
`COBREXA.jl` is quite straightforward:

1. Import the `Distributed` package and add worker processes, e.g. using
   `addprocs`.
2. Pick an analysis function that can be parallelized (such as `screen`
   or `flux_variability_analysis`) and prepare it to work on your data.
3. Pass the desired set of worker IDs to the function using `workers=` argument,
   in the simplest form using e.g. `screen(...,  workers=workers())`.
4. Worker communication will be managed automatically, and you will get results
   "as usual", just appropriately faster.

Specific documentation is available about [running parallel analysis
locally](2_parallel.md) and [running distributed analysis in HPC clusters](3_slurm.md).

## Functions that support parallelization

The functions that support parallel execution include:

- [`flux_variability_analysis`](@ref) (and the underlying [`constraints_variability`](@ref))
- [`screen`](@ref) and [`screen_optimization_model`](@ref)
- [`gene_knockouts`](@ref)
- [`flux_sample`](@ref) (and the underlying [`sample_constraints`](@ref))
- [`objective_production_envelope`](@ref) (and the underlying [`constraints_objective_envelope`](@ref))

Notably, the screening functions can be reused to run many other kinds of
analyses which, in turn, inherit the parallelizability. This includes a wide
range of use-cases that can thus be parallelized very easily:

- single and multiple gene deletions (and other genetic modifications)
- multiple reaction knockouts
- envelope-like production profiles (e.g., enzyme-constrained growth profiles)
- growth media explorations (such as explorations of metabolite depletion)

## Mitigating parallel inefficiencies

Ideally, the speedup gained by parallel processing should be proportional to
the amount of hardware you add as the workers. You should be aware of factors
that reduce the parallel efficiency, which can be summarized as follows:

- Parallelization within single runs of the linear solver is typically not
  supported (and if it is, it may be inefficient for common problem sizes).
  Normally, you want to parallelize the analyzes that comprise multiple
  independent runs of the solvers.
- Some analysis function, such as [`flux_variability_analysis`](@ref), have
  serial parts that can not be parallelized by default. Usually, you may avoid
  the inefficiency by precomputing the serial analysis parts without involving
  the cluster of the workers.
- Frequent worker communication may vastly reduce the efficiency of parallel
  processing; typically this happens if the time required for individual
  analysis steps is smaller than the network round-trip-time to the worker
  processes. Do not use parallelization for very small tasks.
- Transferring large amounts of data among workers may hamper parallel
  efficiency too. Use a single loaded model data object and apply any required
  small modifications directly on the workers to avoid this kind of
  inefficiency.

!!! note "Cost of the distribution and parallelization overhead"
    Before allocating extra resources into the distributed execution, always
    check that your tasks are properly parallelizable and sufficiently large
    to saturate your computation resources, so that the invested energy is not
    wasted.
    [Amdahl's](https://en.wikipedia.org/wiki/Amdahl's_law) and
    [Gustafson's](https://en.wikipedia.org/wiki/Gustafson%27s_law) laws can
    give you a better overview of the sources and consequences of the
    parallelization inefficiencies and the costs of the resulting overhead.
