tools = params.globals.tools
include {ConcatVCF} from './utility'
workflow wf_jointly_genotype_gvcf{
    take: _sample_ids
    take: _gvcfs
    take: _tbis
    take: _fasta
    take: _fasta_fai
    take: _bed_intervls
    take: _target_bed
    take: _dict
    take: _dbsnp
    take: _dbsnp_index
    // take: _bed_intervls
    main:
    
    CombineGVCFs(
        _gvcfs,
        _tbis,
        _fasta,
        _fasta_fai,
        _dict
    )
    // Create a cartesian product of the cohort gvcf and the interval files
    // for parallelization
    ch_int_cohort_gvcf = CombineGVCFs.out.combine(_bed_intervls)
                        //.dump(tag: 'joint_genotype_parallelisation: ')
    // // GenomicsDBImport(mapped_gvcf_GenotypeGVCFs)             
    // _fasta.dump(tag:'fasta')                 
    // _fasta_fai.dump(tag:'fasta_fai')                 
    // _dict.dump(tag:'dict')                 
    // _dbsnp.dump(tag:'dbsnp')                 
    // _dbsnp_index.dump(tag:'_dbsnp_index')
    GenotypeGVCFs(
        ch_int_cohort_gvcf,
        _fasta,
        _fasta_fai,
        _dict,
        _dbsnp,
        _dbsnp_index
    )
    
    // // GenotypeGVCFs output
    // //tuple val("HaplotypeCaller"),  val(patientSampleIdMap), file(interval_bed), file("vcf"), file("vcf.idx")
    
    // // A Cohort vcf
    // vcf_cohort_concatenate_vcf = GenotypeGVCFs.out.vcf_GenotypeGVCFs
    // .map{caller, li_patient_sample_id_map, interval_bed,  vcf, vcf_idx ->
    //     [vcf]
    // }
    // .collect()
    // // .dump(tag: 'cohor_vcf: ')

    CohortConcatVCF(
        GenotypeGVCFs.out[0].collect(),
        GenotypeGVCFs.out[1].collect(),
        _fasta_fai,
        _target_bed,
    )
    _ch_select_variant_input = CohortConcatVCF.out[0].combine(_sample_ids)
    SelectVariants(
        _ch_select_variant_input,
        _fasta,
        _fasta_fai,
        _dict
    )

    emit:
    vcf_with_index = SelectVariants.out[0]
    // vcfs_without_indexes = ConcatVCF.out.concatenated_vcf_without_index
    // cohort_vcf_with_index = CohortConcatVCF.out.cohort_vcf_with_index
    // cohort_vcf_without_index = CohortConcatVCF.out[1]
} // end of wf_haplotypecaller


// STEP GATK HAPLOTYPECALLER.2

// Convert per sample gvcfs collected into a multi-sample gvcf



process CombineGVCFs {
    label 'container_llab'
    label 'cpus_32'
    publishDir "${params.outdir}/VariantCalling/CombinedGVCF", mode: params.publish_dir_mode
    input:
        // tuple val(interval_name), file(interval_bed), val(list_id_patient), val(list_id_sample), file(gdb)
        file(gvcfs)
        file(tbis)
        file(fasta)
        file(fastaFai)
        file(dict)

    output:
       tuple file('cohort.g.vcf.gz'), file('cohort.g.vcf.gz.tbi')
    //    file('cohort.g.vcf.gz.tbi')
    
    when: 'joint_genotype' in tools

    script:
    vcfs_str = ''
    gvcfs.each{vcfs_str += "--variant ${it} "}

    """
    init.sh
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        CombineGVCFs \
        -R ${fasta} \
        ${vcfs_str} \
        --create-output-variant-index \
        -O cohort.g.vcf.gz
    """
}

process GenotypeGVCFs {
    echo true
    label 'container_llab'
    label 'cpus_16'
    tag {intervalBed.baseName}
    input:
        // tuple val(interval_name), file(interval_bed), val(list_id_patient), val(list_id_sample), file(gdb)
        // tuple val(interval_name), file(interval_bed), val(patientSampleIdMap), file(gdb)
        tuple file(cohort_gvcf), file(tbi), file(intervalBed)
        file(fasta)
        file(fastaFai)
        file(dict)
        file(dbsnp)
        file(dbsnpIndex)

    output:
    file("${out_file_bn}.vcf") 
    file("${out_file_bn}.vcf.idx")

    when: 'joint_genotype' in tools

    script:
    out_file_bn = intervalBed.baseName 
    // Using -L is important for speed and we have to index the interval files also
    """
    init.sh
    echo "cohort_gvcf: ${cohort_gvcf}"
    echo "tbi: ${tbi}"
    echo "intervalBed: ${intervalBed}"
    echo "fasta: ${fasta}"
    echo "fastaFai: ${fastaFai}"
    echo "dict: ${dict}"
    echo "dbsnp: ${dbsnp}"
    echo "dbsnpIndex: ${dbsnpIndex}"
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GenotypeGVCFs \
        -R ${fasta} \
        -L ${intervalBed} \
        -D ${dbsnp} \
        -V ${cohort_gvcf} \
        --create-output-variant-index \
        -O "${out_file_bn}.vcf"
    """
}

process CohortConcatVCF {
    label 'container_llab'
    label 'cpus_8'

    tag {'CohortConcatVCF'}

    publishDir "${params.outdir}/VariantCalling/HC_cohort_vcf", mode: params.publish_dir_mode

    input:
        file(vcfs)
        file(tbis)
        file(fastaFai)
        file(targetBED)

    output:
        tuple file("HC_cohort.vcf.gz"), file("HC_cohort.vcf.gz.tbi"), emit: cohort_vcf_with_index
        file("HC_cohort.vcf.gz")

    when: ('haplotypecaller' in tools ||
           'mutect2' in tools || 
           'freebayes' in tools ||
           'joint_genotype' in tools )

    script:
    options = params.target_bed ? "-t ${targetBED}" : ""
    """
    init.sh
    concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o HC_cohort.vcf ${options}
    """
}


process SelectVariants {
    // echo true
    label 'container_llab'
    label 'cpus_4'
    tag {idSample}
    
    publishDir "${params.outdir}/VariantCalling/${idSample}/HC_jointly_genotyped_vcf", mode: params.publish_dir_mode
    input:
        // tuple val(caller), val(id_patient), val(id_sample), val(interval_name), file(interval_bed), file (vcf), file (vcf_idx)
        tuple file (cohort_vcf), file (tbi), idSample
        file(fasta)
        file(fastaFai)
        file(dict)

    output:
    // tuple val("HaplotypeCaller_Jointly_Genotyped"), id_patient, id_sample, file("${interval_bed.baseName}_${id_sample}.vcf"), emit: vcf_SelectVariants
    
    tuple val("HaplotypeCaller_Jointly_Genotyped"), val('patient id placeholder') ,idSample, file("${idSample}.vcf.gz"), file("${idSample}.vcf.gz.tbi")
    tuple file("${idSample}.vcf.gz"), file("${idSample}.vcf.gz.tbi")
    when: ('joint_genotype' in tools )

    script:
    // samples_str = ''
    // sample_ids.each{samples_str+= "${it} "}
    // Using -L is important for speed and we have to index the interval files also
    """
    init.sh
    gatk --java-options -Xmx${task.memory.toGiga()}g \
            SelectVariants \
            -R ${fasta} \
            -V ${cohort_vcf} \
            -O ${idSample}.vcf.gz \
            -sn ${idSample} \
            --create-output-variant-index
    
    """
}






