

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
        ro_input=lambda wildcards, input: utils.as_ro(config, input),
    shell:
        r"""
        tar --create --xz \
            --file {output} \
            --transform 's|.*/\({params.exp}-.*-tier_pdf\..*\)|{params.exp}-pdfs/\1|g' \
            {params.ro_input}
        """


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
        stp_files_regex=utils.as_ro(
            config, patterns.output_simjob_regex(config, tier="stp")
        ),
        ro_evt_files=lambda wildcards, input: utils.as_ro(config, input.evt_files),
    log:
        patterns.log_pdffile_path(config, proctime),
    benchmark:
        patterns.benchmark_pdffile_path(config)
    shell:
        execenv_pyexe(config, "build-pdf") + "--log {log} "
        "-c {input.config_file} "
        "-m $_/inputs "
        "-r {params.stp_files_regex} "
        "-o {output} "
        "-- {input.evt_files}"
