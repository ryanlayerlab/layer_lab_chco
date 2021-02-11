include {reduceVCF} from './utility'
// STEP VCF.QC
tools =  params.globals.tools
skip_qc =  params.globals.skip_qc
workflow wf_vcf_stats{
    // deepvariant output
    //tuple val('DeepVariant'), idSample, file("${idSample}.vcf.gz")
    take: _deepvariant_vcfs
    take: _haplotypecaller_vcfs
    main:
        combined_vcfs = _deepvariant_vcfs
                        .mix(_haplotypecaller_vcfs)
                        // .dump(tag: 'vcfs_stats:')
        BcftoolsStats(combined_vcfs)
        Vcftools(combined_vcfs)
    emit:
        bcfootls_stats = BcftoolsStats.out
        vcfootls_stats = Vcftools.out
}

process BcftoolsStats {
    label 'container_llab'
    label 'cpus_8'

    tag {"${variantCaller} - ${vcf}"}

    publishDir "${params.outdir}/Reports/${idSample}/BCFToolsStats/${variantCaller}", mode: params.publish_dir_mode
    publishDir "${params.outdir}/QC/${idSample}/BCFTools/${variantCaller}", mode: params.publish_dir_mode

    input:
        // tuple variantCaller, idSample, file(vcf)
        tuple variantCaller, idPatient, idSample, file(vcf) , file(vcf_tbi)

    output:
        // file ("*.bcf.tools.stats.out"), emit: bcftools_reports
        file ("*.bcf.tools.stats.out")

    when: !('bcftools' in skip_qc)

    script:
    """
    init.sh
    bcftools stats ${vcf} > ${reduceVCF(vcf.fileName)}.bcf.tools.stats.out
    """
}


process Vcftools {
    label 'container_llab'
    label 'cpus_8'

    tag {"${variantCaller} - ${vcf}"}

    publishDir "${params.outdir}/Reports/${idSample}/VCFTools/${variantCaller}", mode: params.publish_dir_mode
    publishDir "${params.outdir}/QC/${idSample}/VCFTools/${variantCaller}", mode: params.publish_dir_mode

    input:
        tuple variantCaller, idPatient, idSample, file(vcf) , file(vcf_tbi)

    output:
        file ("${reduceVCF(vcf.fileName)}.*")

    when: !('vcftools' in skip_qc)

    script:
    """
    init.sh

    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-count \
    --out ${reduceVCF(vcf.fileName)}

    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-qual \
    --out ${reduceVCF(vcf.fileName)}

    vcftools \
    --gzvcf ${vcf} \
    --FILTER-summary \
    --out ${reduceVCF(vcf.fileName)}
    """
}

