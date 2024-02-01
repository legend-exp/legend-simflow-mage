# legend-simflow

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/legend-exp/legend-simflow?logo=git)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![GitHub issues](https://img.shields.io/github/issues/legend-exp/legend-simflow?logo=github)
![GitHub pull requests](https://img.shields.io/github/issues-pr/legend-exp/legend-simflow?logo=github)
![License](https://img.shields.io/github/license/legend-exp/legend-simflow)

End-to-end Snakemake workflow to run Monte Carlo simulations of signal and background
signatures in the LEGEND experiment and produce probability-density functions (pdfs).
Configuration metadata (e.g. rules for generating simulation macros or
post-processing settings) is stored at
[legend-simflow-config](https://github.com/legend-exp/legend-simflow-config).

## Key concepts

- Simulations are labeled by an unique identifier (e.g. `l200a-hpge-bulk-2vbb`), often
  referred as `simid` (simID). The identifier is defined in
  [legend-simflow-config](https://github.com/legend-exp/legend-simflow-config)
  through `simconfig.json` files in tier directories `raw` and `ver`.
- Each simulation is defined by a template macro (also stored as metadata) and by a set of rules
  (in `simconfig.json`) needed to generate the actual macros (template variable substitutions,
  number of primaries, number of jobs, etc).
- The production is organized in tiers. The state of a simulation in a certain tier is labeled
  as `<tier>.<simid>`. Snakemake understands this syntax.
- The generated pdfs refer to a user-defined selection of LEGEND data taking runs. Such a list
  of runs is specified through the configuration file.
- The production can be restricted to a subset of simulations by passing a list of identifiers
  to Snakemake.

### Workflow steps (tiers)

<p align="center">
  <img src=".github/dag-example.png" alt="DAG example">
</p>

*Figure: representation of an example workflow for simID `l200a-hpge-bulk-2vbb`, 5 parallel jobs
and runIDs `l200-p04-r00{2,3}-phy`. The first word in each box is the name of the corresponding
Snakemake rule.*

1. Macro files are generated and writted to disk according to rules
   defined in the metadata. These macros are in one-to-one correspondence
   with simulation jobs in the next tiers
1. Tier `ver` building: run simulations that generate Monte Carlo
   event vertices needed to some simulations in the next tier. Simulations that do not need
   a special event vertices will directly start from tier `raw`.
1. Tier `raw` building: run full event simulations.
1. Tier `hit` building: run the first (hit-oriented) step of simulation post-processing. Here, "hit"
   stands for Geant4 "step". In this tier, step-wise operations like optical map application
   or step clustering are typically applied.
1. Tier `evt` building: multiple operations are performed in order to build actual events and
   incorporate information about the data taking runs for which the user wants to build pdfs:
   - Partition the `hit` event statistics into fractions corresponding to the actual total livetime
     fraction spanned by each selected run. This information is extracted from
     [`legend-metadata/dataprod/runinfo.json`](https://github.com/legend-exp/legend-metadata/blob/main/dataprod/runinfo.json)
   - Apply HPGe energy resolution functions for each selected run found in the data production
     auto-generated metadata
   - Apply HPGe status flags (available in
     [`legend-metadata/hardware/config`](https://github.com/legend-exp/legend-metadata/blob/main/hardware/config))
1. Tier `pdf` building: summarize `evt`-tier output into histograms (the pdfs).

## Setup

1. Set up [legend-prodenv](https://github.com/legend-exp/legend-prodenv)
   to collect software dependencies, instantiate and manage multiple production
   environments. Snakemake can be installed by following instructions in
   [`legend-prodenv/README.md`](https://github.com/legend-exp/legend-prodenv/blob/main/README.md).

1. Instantiate a new production cycle:
   ```console
   > git clone git@github.com:legend-exp/legend-prodenv
   > cd legend-prodenv && source setup.sh
   > simprod-init-cycle <path-to-cycle-directory>
   ```

1. Customize `<path-to-cycle-directory>/config.json` (see next section).

## The configuration file

The `config.json` file in the production directory allows to customize the workflow
in great detail. Here's a basic description of its fields:

- `experiment`: labels the experiment to be simulated. The same name is used in
  the metadata to label the corresponding configuration files.
- `simlist`: list of simulation identifiers (see below) to be processed by Snakemake.
  Can be a list of strings or a path to a text file. If `*` or `all`, will process all simulations
  defined in the metadata.
- `runlist`: list of LEGEND data taking runs to build pdfs for, in the standard format
  `<experiment>-<period>-<run>-<type>` (e.g. `l200-p03-r000-phy`)
- `benchmark`: section used to configure a benchmarking run:
  - `enabled`: boolean flag to enable/disable the feature
  - `n_primaries`: number of primary events to be simulated in the lower tiers `ver` and `raw`.
- `paths`: customize paths to input or output files.
- `filetypes`: customize the file extensions to be used for input and output files
  generated ny Snakemake.
- `runcmd`: defines for each tier the command to be executed in order to
  produce output files (with Snakemake wildcards, see the corresponding rule for more
  details). The `scripts/MaGe.sh` wrapper, which detects problems in the
  simulation output and exits with a proper code, should be used (and kept up-to-date)
  instead of the `MaGe` command. This allows Snakemake to better detect job failures.
- `execenv`: defines the software environment (container) where all jobs
  should be executed (see below).

> [!TIP]
> all these configuration parameters can be overridden at runtime through Snakemake's `--config` option.

## Production

Run a production by using one of the provided site-specific profiles (recommended):

```console
> cd <path-to-cycle-directory>
> snakemake --profile workflow/profiles/<profile-name>
```

If no system-specific profiles are provided, the `--profile` option can be omitted.
Snakemake will use the `default` profile.

```console
> snakemake
```

The `--config` command line option is very useful to override configuration values.
It can be used, for example, to restrict the production to a subset of simulations:

```console
> snakemake --config simlist="mylist.txt" [...]
```

where `mylist.txt` is a text file in the format:
```
raw.l200a-fibers-Ra224-to-Pb208
hit.l200a-hpge-bulk-2vbb
...
```

One can even just directly pass a comma-separated list:

```console
> snakemake --config simlist="raw.l200a-fibers-Ra224-to-Pb208,hit.l200a-hpge-bulk-2vbb"
```

Once the production is over, the `print_stats` rule can be used to display a table with runtime statistics:

```console
> snakemake -q all print_stats
                                                          wall time [s]         wall time [s]
simid                                                      (cumulative)   jobs      (per job)  primaries
-----                                                     -------------   ----  -------------  ---------
hit.l200a-fibers-Ra224-to-Pb208                                83:20:00    100        0:50:00   1.00E+08
raw.l200a-fibers-Ra224-to-Pb208                                58:20:35    100        0:35:00   1.00E+08
raw.l200a-fibers-Rn222-to-Po214                                33:20:00    100        0:20:00   1.00E+08
...                                                                 ...    ...            ...        ...
```

The `inspect_simjob_logs` rule allows to inspect the simulation log files in search for warnings:
```console
> snakemake inspect_simjob_logs
```
This can generate a lot of output, consider piping it to a file.

Find some useful Snakemake command-line options at the bottom of this page.

> [!IMPORTANT]
> Geant4 macro files are marked as
> [`ancient`](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#ignoring-timestamps)
> as a workaround to the fact that they might have been re-generated (i.e. they have a more recent
> creation time) but with the same content as before (i.e. there is no need
> to re-run the simulations). Since Snakemake removes output files before re-generating,
> there is no simple way to avoid overwriting files with unchanged content.
> The workaround lets user increase simulation statistics by running new jobs (identical to
> those already on disk). However, **If the macro content is updated in the metadata, users will need
> to manually remove the output simulation files or force execution.**

### Benchmarking runs

This workflow implements the possibility to run special "benchmarking" runs in
order to evaluate the speed of simulations, for tuning the number of events to
simulate for each simulation run.

1. Create a new, dedicated production cycle (see above)
2. Enable benchmarking in the configuration file and customize further settings
3. Start the production as usual.

Snakemake will spawn a single job (with the number of primary events specified in the
configuration file) for each simulation. Once the production is over, the results can
be summarized via the `print_benchmark_stats` rule:

```console
> snakemake -q all print_benchmark_stats
simid                                             CPU time [ms/ev]  evts / 1h  jobs (1h) / 10^8 evts
-----                                             ----------------  ---------  ---------------------
raw.l200a-birds-nest-K40                                (13s) 2.79    1288475                     77
raw.l200a-birds-nest-Ra224-to-Pb208                   (191s) 38.33      93916                   1064
raw.l200a-fiber-support-copper-Co60                   (223s) 44.69      80558                   1241
...                                                            ...        ...                    ...
```

> [!NOTE]
> The CPU time is a good measure of the actual simulation time, since other tasks (e.g.
> application loading) are typically not CPU intensive.

## NERSC-specific instructions

### Setup

As an alternative to installing Snakemake through legend-prodenv's tools,
[NERSC's Mamba can be
used](https://docs.nersc.gov/jobs/workflow/snakemake/#building-an-environment-containing-snakemake).

```console
> module load python
> export SWPREFIX="/global/common/software/m2676/<your user>"
> mamba create --prefix $SWPREFIX/.conda/snakemake -c conda-forge -c bioconda snakemake uproot panoptes-ui root
> mamba activate $SWPREFIX/.conda/snakemake
```

> [!IMPORTANT]
> To make proper use of LEGEND's shifter containers, special permissions must be
> set on the `input/simprod/config/MaGe` directory **and all its parents** (see
> [docs](https://docs.nersc.gov/development/shifter/faq-troubleshooting/#invalid-volume-map)):
> ```console
> > setfacl -m u:nobody:x <path-to-cycle-directory>/input/simprod/MaGe
> > setfacl -m u:nobody:x <path-to-cycle-directory>/input/simprod
> > setfacl -m u:nobody:x <path-to-cycle-directory>/input
> > setfacl -m u:nobody:x <path-to-cycle-directory>
> ...and further back if needed...
> ```
> This is not needed if hosting the production below `$PSCRATCH` or if using
> [Podman-HPC](https://docs.nersc.gov/development/podman-hpc/overview).

### Production

Start the production on the interactive node (the default profile works fine):
```
snakemake
```

Start the production on the batch nodes (via SLURM):
```
snakemake --profile workflow/profiles/nersc-batch
```

> [!WARNING]
> This profile does not work as expected at the moment, see https://github.com/legend-exp/legend-simflow/issues/8.
> [This temporary script](https://github.com/legend-exp/legend-simflow/blob/main/profiles/nersc-batch/nersc-submit.sh)
> can be used instead. Note that the maximum runtime at NERSC is 12 hours, so jobs might be killed.
>
> `nersc-submit.sh` usage:
>
> Send a SLURM job for each simulation ID up to the `pdf` tier (in parallel).
> ```console
> > cd <production dir>
> > ./workflow/profiles/nersc-batch/nersc-submit.sh parallel
> ```
> Dry run, to check what would be submitted (recommended):
> ```console
> > DRY_RUN=1 ./workflow/profiles/nersc-batch/nersc-submit.sh
> ```
> Send a Snakemake execution as a single job:
> ```console
> > ./workflow/profiles/nersc-batch/nersc-submit.sh [SNAKEMAKE ARGS]
> ```

## Useful Snakemake CLI options

```
usage: snakemake [OPTIONS] -- [TARGET ...]

  --dry-run, --dryrun, -n
                        Do not execute anything, and display what would be done. If you have a very large workflow,
                        use --dry-run --quiet to just print a summary of the DAG of jobs.
  --profile PROFILE     Name of profile to use for configuring Snakemake.
  --cores [N], -c [N]   Use at most N CPU cores/jobs in parallel. If N is omitted or 'all', the limit is set to the
                        number of available CPU cores.
  --config [KEY=VALUE ...], -C [KEY=VALUE ...]
                        Set or overwrite values in the workflow config object.
  --configfile FILE [FILE ...], --configfiles FILE [FILE ...]
                        Specify or overwrite the config file of the workflow.
  --touch, -t           Touch output files (mark them up to date without really changing them) instead of running
                        their commands.
  --forceall, -F        Force the execution of the selected (or the first) rule and all rules it is dependent on
                        regardless of already created output.
  --list, -l            Show available rules in given Snakefile. (default: False)
  --list-target-rules, --lt
                        Show available target rules in given Snakefile.
  --summary, -S         Print a summary of all files created by the workflow.
  --printshellcmds, -p  Print out the shell commands that will be executed. (default: False)
  --quiet [{progress,rules,all} ...], -q [{progress,rules,all} ...]
                        Do not output certain information. If used without arguments, do not output any progress or
                        rule information. Defining 'all' results in no information being printed at all.
```

## Running jobs in the LEGEND software container

### with [container-env](https://github.com/oschulz/container-env)

In `config.json`:
```json
"execenv": [
    "MESHFILESPATH=$_/inputs/simprod/MaGe/data/legendgeometry/stl_files",
    "MAGERESULTS=$_/inputs/simprod/MaGe/data/legendgeometry",
    "cenv", "legendsw"
]
```
where the `legendsw` environment has been previously created with
```
cenv --create legendsw <path-to-container>
```
and the environment variables are needed for MaGe to work.

### with NERSC Shifter

In `config.json`:
```json
"execenv": [
    "shifter",
    "--image legendexp/legend-software:latest",
    "--volume $_/inputs/simprod/MaGe:/private",
    "--env MESHFILESPATH=/private/data/legendgeometry/stl_files",
    "--env MAGERESULTS=/private/data/legendgeometry"
]
```

### With NERSC Podman-HPC

```js
"execenv": [
  "podman-hpc run",
  "--volume $_/inputs/simprod/MaGe:/private", // mount private MaGe resources
  "--env MESHFILESPATH=/private/data/legendgeometry/stl_files",
  "--env MAGERESULTS=/private/data/legendgeometry",
  "--volume $_:$_", // make production folder available in the container
  "--volume $PSCRATCH:$PSCRATCH", // make scratch area visible too
  "--workdir $$PWD", // podman-hpc does not automatically cd into cwd, unfortunately. NOTE: double $$
  "docker.io/legendexp/legend-software:latest",
  "--"
]
```
