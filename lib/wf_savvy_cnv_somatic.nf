tools =  params.globals.tools
workflow wf_savvy_cnv_somatic{
    take: _bam_recal
    take: _savvy_controls_dir
    take: _savvy_sizes
    main:
     /* SavvyCNV Somatic Copy Number related calls */
    /* Starting point is duplicated marked bams from MarkDuplicates.out.marked_bams with the following structure */
    /* MarkDuplicates.out.marked_bams => [idPatient, idSample, md.bam, md.bam.bai]*/
    SavvyCNVCoverageSummary(_bam_recal)
    SavvyCNV(SavvyCNVCoverageSummary.out.collect(),
             _savvy_controls_dir)

    SavvyCNV_params(SavvyCNVCoverageSummary.out.collect(), _savvy_controls_dir, _savvy_sizes)

    emit:
    savvy_output = SavvyCNV.out
    savvy_param_output = SavvyCNV_params.out
} // end of wf_germline_cnv


// /* SavvyCNV Related processes */

process SavvyCNVCoverageSummary {
    label 'container_llab'
    label 'cpus_16'
    tag "${idSample}"
    // publishDir "${params.outdir}/VariantCalling/${idSample}/SavvycnvCoverageSummary", mode: params.publish_dir_mode
    
    input:
         tuple idPatient, idSample, file(bam), file(bai)
    
    output:
    file("${idSample}.coverageBinner")

    when: 'savvy_cnv_somatic' in tools

    script:
    
    """
    init.sh
    java -Xmx1g CoverageBinner ${bam} > ${idSample}.coverageBinner
    """
}

process SavvyCNV {
    label 'container_llab'
    label 'cpus_16'
    // tag "${idSample}"
    publishDir "${params.outdir}/VariantCalling/", mode: params.publish_dir_mode
    
    input:
        file("*")
        file(savvy_controls_dir)
     
    output:
        file("SavvycnvResults")
    
    when: 'savvy_cnv_somatic' in tools
    
    script:
    //chunk_size = 200000
    chunk_size = 250
    control_options = params.savvy_controls_dir ? " -control ${savvy_controls_dir}/*.coverageBinner ": ""
    """
    init.sh
    mkdir -p SavvycnvResults/SavvycnvCoverageSummary
    mkdir  SavvycnvResults/pdfs
    java -Xmx30g SavvyCNV -a -d ${chunk_size} -case *.coverageBinner ${control_options} > cnv_list.bed 2>log_messages.txt
    cp *.coverageBinner SavvycnvResults/SavvycnvCoverageSummary/
    cp *.cnvs.pdf SavvycnvResults/pdfs
    cp cnv_list.bed log_messages.txt SavvycnvResults/
    """
}

process SavvyCNV_params {
    label 'container_llab'
    label 'cpus_16'
    // tag "${idSample}"
    publishDir "${params.outdir}/VariantCalling/SavvySize${savvy_size}", mode: params.publish_dir_mode

    input:
        file("*")
        file(savvy_controls_dir)
        val savvy_size

    output:
        tuple file("SavvycnvResults${savvy_size}"), file("cnv_list_${savvy_size}.bed")

    when: 'savvy_cnv_somatic' in tools

    script:
    //chunk_size = 200000
    chunk_size = 250
    control_options = params.savvy_controls_dir ? " -control ${savvy_controls_dir}/*.coverageBinner ": ""
    """
    init.sh
    mkdir -p SavvycnvResults${savvy_size}/SavvycnvCoverageSummary
    mkdir  SavvycnvResults${savvy_size}/pdfs
    java -Xmx30g SavvyCNV -a -d ${savvy_size} -case *.coverageBinner ${control_options} > cnv_list_${savvy_size}.bed 2>log_messages.txt
    cp *.coverageBinner SavvycnvResults${savvy_size}/SavvycnvCoverageSummary/
    cp *.cnvs.pdf SavvycnvResults${savvy_size}/pdfs
    cp cnv_list_${savvy_size}.bed log_messages.txt SavvycnvResults${savvy_size}/
    """
}
