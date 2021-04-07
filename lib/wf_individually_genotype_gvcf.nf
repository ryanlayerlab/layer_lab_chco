tools = params.globals.tools
include {ConcatVCF} from './utility'
workflow wf_individually_genotype_gvcf{
    // input structure
    // tuple val("HaplotypeCallerGVCF"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.g.vcf")
    take: _gvcf_HC
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _dbsnp
    take: _dbsnp_index
    take: _target_bed
    main:
        IndividuallyGentoypeGVCF(
            _gvcf_HC,
            _fasta,
            _fasta_fai,
            _dict,
            _dbsnp,
            _dbsnp_index
        )
        
        _vcf_concatVCF = 
            IndividuallyGentoypeGVCF.out.vcf_HaplotypeCaller
            .groupTuple(by: [0,1])
            .map{ idPatient, idSample, vcfs -> 
                ['HaplotypeCaller_Individually_Genotyped', idPatient, idSample, vcfs]
            }
            // .dump(tag: "vcf_ConcatVcf")

        ConcatVCF(
                _vcf_concatVCF,
                _fasta_fai,
                _target_bed,
                'HC', // prefix for output files
                'vcf', // extension for the output files
                'HC_individually_genotyped_vcf' // output directory name
                )
        // GvcfToVcf(
        //         ConcatVCF.out.concatenated_vcf_without_index
        //         )
    emit:
        sample_vcf_HC = ConcatVCF.out.concatenated_vcf_with_index
        //sample_vcf_HC = IndividuallyGentoypeGVCF.out.vcf_HaplotypeCaller
        // sample_vcf_HC = GvcfToVcf.out.vcf_HaplotypeCaller
}

process IndividuallyGentoypeGVCF{
    label 'container_llab'
    label 'cpus_8'
    tag {idSample + "-" + gvcf.baseName}
    // tag {idSample} 
    // publishDir "${params.outdir}/VariantCalling/${idSample}/HC_individually_genotyped_vcf", mode: params.publish_dir_mode
    input:
        tuple idPatient, idSample, file(intervalsBed), file(gvcf)
        file(fasta)
        file(fastaFai)
        file(dict)
        file(dbsnp)
        file(dbsnpIndex)
    output:
        //tuple val('HaplotypeCaller_Individually_Genotyped'), idPatient, idSample, file("${gvcf.simpleName}.vcf"), file("${gvcf.simpleName}.vcf.tbi"), emit: vcf_HaplotypeCaller
        tuple idPatient, idSample, file(out_file), emit: vcf_HaplotypeCaller

    when: 'haplotypecaller' in tools

    script:
    // fn=gvcf.fileName
    // prefix=fn.minus(".g.vcf")
    // out_file="${gvcf.fileName}.vcf"
    prefix="${gvcf.fileName}" - ".g.vcf"
    // We'll first decompress the gvcf
    // in_file= "${gvcf.fileName}" - ".gz"
    out_file="${prefix}.vcf"
    """
    init.sh
    bgzip  ${gvcf}
    tabix  ${gvcf}.gz
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GenotypeGVCFs \
        -R ${fasta} \
        -L ${intervalsBed} \
        -D ${dbsnp} \
        -V ${gvcf}.gz \
        -O "${out_file}"
    """
}
