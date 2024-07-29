
# Local parallel processing

To run an analysis in parallel, we first need to load the `Distributed`
package and add a few worker processes. For example, we may start 5 local
processes (that may utilize 5 CPUs) as follows

```julia
using Distributed
addprocs(5)
```

!!! note "`Distributed.jl` installation"
    `Distributed.jl` usually comes pre-installed with Julia distribution, but
    one may still need to "enable" it by typing `] add Distributed`.

To check that the workers are really there, use `workers()`. In this
case, it should return a vector of _worker IDs_, very likely equal to
`[2,3,4,5,6]`.

Each of the processes contains a self-sufficient image of Julia that can act
independently; in turn the additional processes also consume some memory. Each
process with loaded `COBREXA.jl` and a simple solver such as HiGHS may consume
around 500MB of RAM, which should be taken into account when planning the
analysis scale.

!!! warning "Using Julia environments with Distributed"
    In certain conditions, the Distributed package does not properly forward
    the project configuration to the workers, resulting to package version
    mismatches and other problems. For pipelines that run in custom project
    folders, use the following form of `addprocs` instead:
    ```julia
    addprocs(5, exeflags=`--project=$(Base.active_project())`)
    ```

Packages (COBREXA and the selected solver) must be loaded at all processes,
which may ensured using the "everywhere" macro (from `Distributed` package):
```julia
@everywhere using COBREXA, HiGHS
```

Utilizing the prepared worker processes is then straightforward: We pass the
list of workers to the selected analysis function using the `workers` keyword
argument, and the parallel processing is orchestrated automatically:

```julia
model = load_model("e_coli_core.xml")
result = flux_variability_analysis(
    model,
    optimizer = HiGHS.Optimizer,
    workers = workers()
)
```
