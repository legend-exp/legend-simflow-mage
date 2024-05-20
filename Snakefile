# Copyright (C) 2023 Luigi Pertoldi <gipert@pm.me>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from pathlib import Path

from scripts.utils import utils, patterns, aggregate, tier_evt

if not config:
    raise RuntimeError("you must set a config file with --configfile")

utils.subst_vars_in_snakemake_config(workflow, config)

config.setdefault("benchmark", {"enabled": False})


wildcard_constraints:
    tier=r"\w+",
    simid=r"[-\w]+",
    jobid=r"\d+",
    runid=r"[-\w]+",


rule gen_all_macros:
    """Aggregate and produce all the macro files."""
    input:
        aggregate.gen_list_of_all_macros(config, tier="ver"),
        aggregate.gen_list_of_all_macros(config, tier="raw"),


rule gen_all_tier_raw:
    """Aggregate and produce all the raw tier files."""
    input:
        aggregate.gen_list_of_all_plots_outputs(config, tier="raw"),
        aggregate.gen_list_of_all_simid_outputs(config, tier="raw"),


rule gen_all_tier_hit:
    """Aggregate and produce all the hit tier files."""
    input:
        aggregate.gen_list_of_all_simid_outputs(config, tier="hit"),


rule gen_all_tier_evt:
    """Aggregate and produce all the evt tier files."""
    input:
        aggregate.gen_list_of_all_tier_evt_outputs(config),


rule gen_all_tier_pdf:
    """Aggregate and produce all the pdf tier files."""
    input:
        aggregate.gen_list_of_all_tier_pdf_outputs(config),


rule gen_pdf_release:
    """Generates a tarball with all the pdf files."""
    message:
        "Generating pdf release"
    input:
        aggregate.gen_list_of_all_tier_pdf_outputs(config),
    output:
        Path(config["paths"]["pdf_releases"]) / (config["experiment"] + "-pdfs.tar.xz"),
    params:
        exp=config["experiment"],
    shell:
        r"""
        tar --create --xz \
            --file {output} \
            --transform 's|.*/\({params.exp}-.*-tier_pdf\..*\)|{params.exp}-pdfs/\1|g' \
            {input}
        """


def gen_target_all():
    if config.get("simlist", "*") in ("all", "*"):
        return (
            rules.gen_pdf_release.output,
            [
                aggregate.gen_list_of_all_plots_outputs(config, tier=t)
                for t in ("ver", "raw")
            ],
        )
    else:
        return aggregate.process_simlist(config)


rule all:
    default_target: True
    input:
        gen_target_all(),


# since the number of generated macros for the 'output' field
# must be deduced at runtime from the JSON configuration, we need here to
# generate a separate rule for each 'simid'
simconfigs = aggregate.collect_simconfigs(config, ["ver", "raw"])

for tier, simid, n_macros in simconfigs:

    rule:
        f"""Generates all needed simulation macros ({n_macros}) for {simid} in tier {tier}. No wildcards are used."""
        localrule: True
        input:
            **patterns.macro_gen_inputs(config, tier, simid),
        output:
            patterns.input_simid_filenames(config, n_macros, tier=tier, simid=simid),
        params:
            tier=tier,
            simid=simid,
        threads: 1
        message:
            f"Generating macros for {tier}.{simid}"
        script:
            "scripts/generate_macros.py"

    utils.set_last_rule_name(workflow, f"gen_macros_{simid}-tier_{tier}")


rule build_tier_ver:
    """Run a single simulation job for the ver tier.
    Uses wildcards `simid` and `jobid`.

    Warning
    -------
    The macro file is marked as "ancient" as a workaround to the fact that
    it might have been re-generated (i.e. it effectively has a more recent
    creation time) but with the same content as before (i.e. there is no need
    to re-run the simulation). If the macro content is updated, users will need
    to manually remove the output simulation files or force execution.
    """
    message:
        "Producing output file for job ver.{wildcards.simid}.{wildcards.jobid}"
    input:
        macro=ancient(patterns.input_simjob_filename(config, tier="ver")),
    output:
        protected(patterns.output_simjob_filename(config, tier="ver")),
    log:
        patterns.log_file_path(config, tier="ver"),
    benchmark:
        patterns.benchmark_file_path(config, tier="ver")
    shell:
        patterns.run_command(config, "ver")


rule build_tier_raw:
    """Run a single simulation job for the raw tier.
    Uses wildcards `simid` and `jobid`.

    Warning
    -------
    The macro file is marked as "ancient" as a workaround to the fact that
    it might have been re-generated (i.e. it effectively has a more recent
    creation time) but with the same content as before (i.e. there is no need
    to re-run the simulation). If the macro content is updated, users will need
    to manually remove the output simulation files or force execution.
    """
    message:
        "Producing output file for job raw.{wildcards.simid}.{wildcards.jobid}"
    input:
        macro=ancient(patterns.input_simjob_filename(config, tier="raw")),
        verfile=lambda wildcards: patterns.smk_ver_filename_for_raw(config, wildcards),
    output:
        protected(patterns.output_simjob_filename(config, tier="raw")),
    log:
        patterns.log_file_path(config, tier="raw"),
    benchmark:
        patterns.benchmark_file_path(config, tier="raw")
    shell:
        patterns.run_command(config, "raw")


rule build_tier_hit:
    """Produces a hit tier file starting from a single raw tier file."""
    message:
        "Producing output file for job hit.{wildcards.simid}.{wildcards.jobid}"
    input:
        raw_file=rules.build_tier_raw.output,
        optmap_lar=config["paths"]["optical_maps"]["lar"],
        optmap_pen=config["paths"]["optical_maps"]["pen"],
        optmap_fiber=config["paths"]["optical_maps"]["fiber"],
    output:
        patterns.output_simjob_filename(config, tier="hit"),
    log:
        patterns.log_file_path(config, tier="hit"),
    benchmark:
        patterns.benchmark_file_path(config, tier="hit")
    shell:
        patterns.run_command(config, "hit")


rule make_tier_evt_config_file:
    """Generates configuration files for `build_tier_evt` based on metadata.
    Uses wildcard `runid`."""
    localrule: True
    input:
        # FIXME: need to list actual files, not the directory
        config["paths"]["metadata"],
    output:
        Path(config["paths"]["genconfig"]) / "{runid}-build_evt.json",
    script:
        "scripts/make_tier_evt_config_file.py"


rule make_run_partition_file:
    """Computes and stores on disk rules for partitioning the simulated event
    statistics according to data taking runs. Uses wildcard `simid`."""
    localrule: True
    input:
        hit_files=lambda wildcards: aggregate.gen_list_of_simid_outputs(
            config, tier="hit", simid=wildcards.simid
        ),
        runinfo=Path(config["paths"]["metadata"]) / "dataprod" / "runinfo.json",
    output:
        Path(config["paths"]["genconfig"]) / "{simid}-run_partition.json",
    script:
        "scripts/make_run_partition_file.py"


rule build_tier_evt:
    """Produces an evt tier file."""
    message:
        "Producing output file for job evt.{wildcards.simid}.{wildcards.runid}"
    input:
        hit_files=lambda wildcards: aggregate.gen_list_of_simid_outputs(
            config, tier="hit", simid=wildcards.simid
        ),
        config_file=rules.make_tier_evt_config_file.output,
        run_part_file=rules.make_run_partition_file.output,
        hpge_db=Path(config["paths"]["metadata"])
        / "hardware/detectors/germanium/diodes",
    output:
        patterns.output_evt_filename(config),
    params:
        evt_window=lambda wildcards, input: tier_evt.smk_get_evt_window(
            wildcards, input
        ),
        hit_files_regex=patterns.output_simjob_regex(config, tier="hit"),
    log:
        patterns.log_evtfile_path(config),
    benchmark:
        patterns.benchmark_evtfile_path(config)
    shell:
        patterns.run_command(config, "evt")


rule build_tier_pdf:
    """Produces a pdf tier file."""
    message:
        "Producing output file for job pdf.{wildcards.simid}"
    input:
        evt_files=lambda wildcards: aggregate.gen_list_of_tier_evt_outputs(
            config, wildcards.simid
        ),
        config_file=patterns.pdf_config_path(config),
    output:
        patterns.output_pdf_filename(config),
    params:
        raw_files_regex=patterns.output_simjob_regex(config, tier="raw"),
    log:
        patterns.log_pdffile_path(config),
    benchmark:
        patterns.benchmark_pdffile_path(config)
    shell:
        patterns.run_command(config, "pdf")


include: "rules/aux.smk"
