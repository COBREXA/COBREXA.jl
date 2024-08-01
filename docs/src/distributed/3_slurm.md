
# Working in a HPC environment

Many researchers have access to institutional HPC facilities that allow
time-sharing of the capacity of a large computer cluster between many users.
Julia and `COBREXA.jl` work well within this environment, and the COBREXA
analyses usually require only minimal additional customization to be able to
find and utilize the resources available from the HPC.

When executed in a HPC environment, the analysis script must solve several
relatively complex tasks:

- It needs to find out how many resources were allocated for the analysis
- It needs to add the remote workers precisely at the allocated places

Fortunately, the package
[`ClusterManagers.jl`](https://github.com/JuliaParallel/ClusterManagers.jl)
does that for us. For simplicity, here we assume that the HPC is scheduled by
[Slurm](https://slurm.schedmd.com/), but other scheduling environments are
supported in a very similar way.

## Interacting with Slurm

Utilization of the Slurm-provided resources is enabled as follows:
- first, import the `ClusterManagers` package
- find how many processes to spawn from the environment, typically from
  `SLURM_NTASKS` environment variable
- use the function `addprocs_slurm` to precisely connect to the allocated
  computational resources

After adding the Slurm workers, one may continue as if the workers were added
using normal `addprocs` --- typically, we can load the model and (for example) run
the `flux_variability_analysis` as if we would use the [local
workers](2_parallel.md).

The Julia script that does a parallel analysis in a Slurm cluster may look as
follows:

```julia
using COBREXA, Distributed, ClusterManagers, HiGHS

available_workers = parse(Int, ENV["SLURM_NTASKS"])

addprocs_slurm(available_workers)

# ... load models, prepare data, etc. ...

results = flux_variability_analysis(..., workers=workers())

# ... save the results into a file ...
```

!!! tip "What about the other HPC schedulers?"
    `ClusterManagers.jl` supports many other common HPC scheduling systems,
    including LFS, Sun Grid, SGE, PBS, and Scyld, in a way almost identical to
    Slurm. See the
    [package documentation](https://github.com/JuliaParallel/ClusterManagers.jl/blob/master/README.md)
    for details.

!!! warning "Using Julia environments with Distributed"
    Sometimes the project configuration is not forwarded to the workers
    automatically, resulting to package version mismatches and other problems.
    When utilizing custom project folders (by running Julia with `julia
    --project=...`), use the following form of `addprocs_slurm` instead:
    ```julia
    addprocs_slurm(available_workers, exeflags=`--project=$(Base.active_project())`)
    ```

## Wrapping a pipeline script in a Slurm batch job

To be able to submit a script for later processing using the [`sbatch` Slurm
command](https://slurm.schedmd.com/sbatch.html), we need to wrap it in a small
"batch" script that tells Slurm how many resources the process needs.

Assuming we have a Julia computation script written down in `myJob.jl` and
saved on the HPC cluster's access node, the corresponding Slurm batch script
(let's call it `myJob.sbatch`) may look as follows:

```sh
#!/bin/bash -l
#SBATCH -n 100           # the job will use 100 individual worker processes
#SBATCH -c 1             # each worker will sit on a single CPU
#SBATCH -t 30            # the whole job will take less than 30 minutes
#SBATCH -J myJob         # the name of the job (for own reference)

module load lang/Julia   # add Julia to the environment (this may differ on different clusters and installations!)

julia myJob.jl
```

To run the computation, run `sbatch myJob.sbatch` on the cluster access node.
The job will be scheduled and eventually executed. It is possible to watch the
output of commands `sacct` and `squeue` in the meantime, to see the progress.

Remember that it is necessary to explicitly save the result of the Julia script
computation to files, to be able to retrieve them later. Standard outputs of
the jobs are often mangled and/or discarded. If we would still want to collect
the standard output of the Julia script, we might need to change the last line
of the batch script as follows:

```sh
julia myJob.jl > myJob.log
```

...and collect the output from `myJob.log` later. This is convenient especially
if the script prints out various computation details using `@info` and similar
macros.
