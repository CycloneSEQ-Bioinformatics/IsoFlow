rule generateEvents:
    input:
        gtf = config["reference"][config["specie"]]["gtf"]
    output:
        event = "diff_splice/events.ioe",
        isoform = "diff_splice/isoforms.ioi"
    params:
        outdir = "diff_splice"
    log: "logs/generateEvents.log"
    benchmark: "benchmarks/generateEvents.benchmark"
    shell:"""
    # generate local AS events
    python {SUPPA_HOME}/suppa.py generateEvents -i {input.gtf} -o {params.outdir}/events -e SE SS MX RI FL -f ioe > {log} 2>&1
    awk 'FNR==1 && NR!=1 {{while (/^<header>/) getline;}} {{print}}' {params.outdir}/*.ioe > {params.outdir}/events.ioe

    # generate the transcript "events"
    python {SUPPA_HOME}/suppa.py generateEvents -i {input.gtf} -o {params.outdir}/isoforms -f ioi >> {log} 2>&1
    """

rule trs_quant:
    input:
        count_tsvs = expand("count_reads/{sample}/quant.sf", sample=all_samples.keys())
    output:
        tpm_tsv = "diff_splice/iso.tpm",
    log: "logs/trs_quant.log"
    benchmark: "benchmarks/trs_quant.benchmark"
    shell:"""
    python {SUPPA_HOME}/multipleFieldSelection.py -i {input.count_tsvs} -k 1 -f 4 -o {output.tpm_tsv} > {log} 2>&1
    """

rule calcPSI:
    input:
        event = rules.generateEvents.output.event,
        isoform = rules.generateEvents.output.isoform,
        tpm_tsv = rules.trs_quant.output.tpm_tsv,
        gtf = config["reference"][config["specie"]]["gtf"],
        trs_info=rules.build_transcriptome_index.output.trs_info,
    output:
        event = "diff_splice/events.psi",
        isoform = "diff_splice/iso_isoform.psi"
    params:
        outdir = "diff_splice"
    log: "logs/calcPSI.log"
    benchmark: "benchmarks/calcPSI.benchmark"
    shell:"""
    # compute the PSI values of the events
    python {SUPPA_HOME}/suppa.py psiPerEvent -i {input.event} -e {input.tpm_tsv} -o {params.outdir}/events > {log} 2>&1
    {SNAKEDIR}/scripts/add_gene_name.R {input.trs_info} {output.event} rowname
    # compute the PSI values of the isoforms
    python {SUPPA_HOME}/suppa.py psiPerIsoform -g {input.gtf} -e {input.tpm_tsv} -o {params.outdir}/iso >> {log} 2>&1
    {SNAKEDIR}/scripts/add_gene_name.R {input.trs_info} {output.isoform} rowname
    """

rule diffSplice:
    input:
        tpm_tsv = rules.trs_quant.output.tpm_tsv,
        event_psi = rules.calcPSI.output.event,
        isoform_psi = rules.calcPSI.output.isoform,
        event = rules.generateEvents.output.event,
        isoform = rules.generateEvents.output.isoform,
        trs_info=rules.build_transcriptome_index.output.trs_info,
    output:
        event_dpsi = "diff_splice/suppa_diffSplice_event.dpsi",
        isoform_dpsi = "diff_splice/suppa_diffSplice_iso.dpsi"
    params:
        ctrl_samples = ','.join(manifest.loc[manifest['condition']=='control','alias'].tolist()),
        treat_samples = ','.join(manifest.loc[manifest['condition']!='control','alias'].tolist()),
        outdir = 'diff_splice'
    log: "logs/diffSplice.log"
    benchmark: "benchmarks/diffSplice.benchmark"
    shell:"""
    # Differential splicing with local events
    ## Split the PSI and TPM files between the 2 conditions:
    Rscript {SUPPA_HOME}/scripts/split_file.R {input.tpm_tsv} {params.ctrl_samples} {params.treat_samples} {params.outdir}/iso.ctrl.tpm {params.outdir}/iso.treat.tpm -i > {log} 2>&1
    Rscript {SUPPA_HOME}/scripts/split_file.R {input.event_psi} {params.ctrl_samples} {params.treat_samples} {params.outdir}/events.ctrl.psi {params.outdir}/events.treat.psi -e >> {log} 2>&1
    ## perform the differential splicing analysis
    python {SUPPA_HOME}/suppa.py diffSplice -m classical -gc -i {input.event} -p {params.outdir}/events.ctrl.psi {params.outdir}/events.treat.psi -e {params.outdir}/iso.ctrl.tpm {params.outdir}/iso.treat.tpm -th 1 -o suppa_diffSplice_event >> {log} 2>&1
    
    # Differential transcript usage
    ## Split the PSI files between 2 conditions:
    Rscript {SUPPA_HOME}/scripts/split_file.R {input.isoform_psi} {params.ctrl_samples} {params.treat_samples} {params.outdir}/iso_isoform.ctrl.psi {params.outdir}/iso_isoform.treat.psi -i >> {log} 2>&1
    ## perform the differential splicing analysis
    python {SUPPA_HOME}/suppa.py diffSplice -m classical -gc -i {input.isoform} -p {params.outdir}/iso_isoform.ctrl.psi {params.outdir}/iso_isoform.treat.psi -e {params.outdir}/iso.ctrl.tpm {params.outdir}/iso.treat.tpm -th 1 -o suppa_diffSplice_iso >> {log} 2>&1
    mv suppa_diffSplice_event* {params.outdir}/
    mv suppa_diffSplice_iso* {params.outdir}/
    {SNAKEDIR}/scripts/add_gene_name.R {input.trs_info} {output.event_dpsi} rowname
    {SNAKEDIR}/scripts/add_gene_name.R {input.trs_info} {output.isoform_dpsi} rowname
    """

rule write_bamlist:
    input:
        event_dpsi = rules.diffSplice.output.event_dpsi,
        bams = expand("genome_alignments/{sample}/aligned.sort.bam", sample=all_samples.keys()),
    output:
        bamlist = "diff_splice/plot_sashimi/input_bams.tsv",
    run:
        manifest['bam'] = manifest['alias'].apply(lambda x: os.path.join(os.getcwd(), "genome_alignments", x, "aligned.sort.bam"))
        manifest.loc[:, ["alias", "bam", "condition"]].to_csv(output.bamlist, sep="\t", header=False, index=False)

rule plot_sashimi:
    input:
        bamlist = rules.write_bamlist.output.bamlist,
        event_dpsi = rules.diffSplice.output.event_dpsi,
        gtf = config["reference"][config["specie"]]["gtf"],
        genome = config["reference"][config["specie"]]["fasta"]
    output:
        plot_log = "diff_splice/plot_sashimi/plot.log",
    params:
        outdir = 'diff_splice/plot_sashimi',
        n_samples = len(all_samples.keys()),
        AS_dPSI = config['AS_dPSI'],
        AS_pval = config['AS_pval'],
    log: "logs/plot_sashimi.log"
    benchmark: "benchmarks/plot_sashimi.benchmark"
    conda: "ggsashimi_env"
    shell:"""
    if [ ! -f {input.genome}.fai ];then
        samtools faidx {input.genome}
    fi
    awk -v AS_dPSI={params.AS_dPSI} -v AS_pval={params.AS_pval} -v OFS='\\t' 'NR>1 && ($2<=-AS_dPSI || $2>=AS_dPSI) && ($3!="nan" && $3<=AS_pval){{split($1,a,/;/);split(a[2],b,/:/);chrom=b[2];split(b[3],c,/-/);start=c[1];strand=b[length(b)];split(b[length(b)-1],c,/-/);end=c[length(c)];print chrom,start,end,strand,$0}}' {input.event_dpsi} | sort -k7,7n | bedtools slop -b 1000 -g {input.genome}.fai -i - > {params.outdir}/suppa_diffSplice_event.sig.dpsi
    if [ "$(wc -l < {params.outdir}/suppa_diffSplice_event.sig.dpsi)" -eq 0 ]; then
        echo "There is no siginificant differential splicing event" > {params.outdir}/plot.log
    else
        while IFS=$'\\t' read -r chrom start end strand event_id event_dPSI event_pval;do
            event_id_new=$(echo "$event_id" | sed 's/[:;]/_/g')
            if [ $strand == "+" ];then
                strand=plus
            else
                strand=minus
            fi
            echo ">$event_id" >> {log}
            {SNAKEDIR}/scripts/ggsashimi.py -b {input.bamlist} -c $chrom:$start-$end -S $strand -g {input.gtf} -C 3 -F png --height 1.5 --ann-height 2 --width $((2*{params.n_samples}+2)) -o {params.outdir}/${{event_id_new}}.${{event_dPSI}}.${{event_pval}}.png >> {log} 2>&1
            if (echo $event_id | grep -q SE) && [ -f {params.outdir}/${{event_id_new}}.${{event_dPSI}}.${{event_pval}}.png ] && [ ! -f {params.outdir}/which_to_show ];then
                echo "${{event_id_new}}.${{event_dPSI}}.${{event_pval}}.png" > {params.outdir}/which_to_show
            fi
            if [ -f {params.outdir}/${{event_id_new}}.${{event_dPSI}}.${{event_pval}}.png ] && [ ! -f {params.outdir}/which_most_sig ];then
                echo "${{event_id_new}}.${{event_dPSI}}.${{event_pval}}.png" > {params.outdir}/which_most_sig
            fi
        done < {params.outdir}/suppa_diffSplice_event.sig.dpsi
        echo "Plot successfully" > {params.outdir}/plot.log
    fi
    """
