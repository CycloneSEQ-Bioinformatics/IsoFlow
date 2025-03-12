rule write_coldata:
    input:
        count_tsv = rules.merge_counts.output.count_tsv
    output:
        coldata = "diff_exp/coldata.tsv"
    run:
        manifest.loc[:, ["sample", "alias", "condition"]].to_csv(output.coldata, sep="\t", index=False)

rule write_de_params:
    input:
        count_tsv = rules.merge_counts.output.count_tsv,
        trs_gtf = config["reference"][config["specie"]]["gtf"]
    output:
        de_params = "diff_exp/de_params.tsv"
    run:
        d = OrderedDict()
        d["ref_annotation"] = [input.trs_gtf]
        d["min_samps_gene_expr"] = [config["min_samps_gene_expr"]]
        d["min_samps_feature_expr"] = [config["min_samps_feature_expr"]]
        d["min_gene_expr"] = [config["min_gene_expr"]]
        d["min_feature_expr"] = [config["min_feature_expr"]]
        df = pd.DataFrame(d)
        df.to_csv(output.de_params, sep="\t", index=False)

rule diff_exp:
    input:
        de_params = rules.write_de_params.output.de_params,
        coldata = rules.write_coldata.output.coldata,
        trs_info = rules.build_transcriptome_index.output.trs_info,
        count_tsv = rules.merge_counts.output.count_tsv,
    output:
        dge_tsv = "diff_exp/results_dge.tsv",
        dge_maplot = "diff_exp/results_dge.maplot.png",
        dge_volcanoplot = "diff_exp/results_dge.volcanoplot.png",
        dte_tsv = "diff_exp/results_dte.tsv",
        dte_maplot = "diff_exp/results_dte.maplot.png",
        dte_volcanoplot = "diff_exp/results_dte.volcanoplot.png",
    params:
        outdir = "diff_exp"
    log:
        "logs/diff_exp.log"
    benchmark:
        "benchmarks/diff_exp.benchmark"
    shell:"""
    {SNAKEDIR}/scripts/de_analysis.R {input.de_params} {input.coldata} {input.trs_info} {input.count_tsv} {params.outdir} >{log} 2>&1
    """
