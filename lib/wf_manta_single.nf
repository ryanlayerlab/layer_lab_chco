
tools = params.globals.tools
status_map = params.globals.status_map
workflow wf_manta_single{
    take: _recal_bam
    take: _target_bed
    take: _fasta
    take: _fasta_fai
    
    main:
        MantaSingle(
            _recal_bam,
            _fasta,
            _fasta_fai,
            _target_bed
        )
    emit:
        output_tuple = MantaSingle.out
} // end of wf_germline_cnv


// STEP MANTA.1 - SINGLE MODE

process MantaSingle {
    label 'container_sarek'
    label 'cpus_32'
    // label 'memory_max'

    tag {idSample}

    publishDir "${params.outdir}/VariantCalling/${idSample}/Manta", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam), file(bai)
        file(fasta)
        file(fastaFai)
        file(targetBED)

    output:
        tuple val("Manta"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi")

    when: 'manta_single' in tools

    script:
    beforeScript=""
    // beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    //options = params.target_bed ? "--exome  --callRegions call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome " : ""
    status = status_map[idPatient, idSample]
    inputbam = status == 0 ? "--bam" : "--tumorBam"
    vcftype = status == 0 ? "diploid" : "tumor"
    """
    ${beforeScript}
    configManta.py \
        ${inputbam} ${bam} \
        --reference ${fasta} \
        ${options} \
        --runDir Manta

    python Manta/runWorkflow.py -m local -j ${task.cpus}

    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${idSample}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${idSample}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${idSample}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/${vcftype}SV.vcf.gz \
        Manta_${idSample}.${vcftype}SV.vcf.gz
    mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
        Manta_${idSample}.${vcftype}SV.vcf.gz.tbi
    """
}
