tools = params.globals.tools
status_map = params.globals.status_map
gender_map = params.globals.gender_map
// tools = params.globals.tools
include {ConcatVCF} from './utility'
workflow wf_haplotypecaller{
    take: _int_bam_recal
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _target_bed
    take: _dbsnp
    take: _dbsnp_index
    // take: _bed_intervls
    main:
        HaplotypeCaller(
            _int_bam_recal,
            _fasta,
            _fasta_fai,
            _dict,
            _dbsnp,
            _dbsnp_index
        )
     // Create individual gvcfs without any genotyping
     gvcf_ConcatVCF = 
            HaplotypeCaller.out.gvcf_GenotypeGVCFs
            .groupTuple(by: [0,1])
            .map{ idPatient, idSample, interval_beds,  gvcfs -> 
                ['HaplotypeCaller_gvcf', idPatient, idSample, gvcfs]
            }
            // .dump(tag: "gvcfs_ConcatVcf")

    ConcatVCF(
                gvcf_ConcatVCF,
                _fasta_fai,
                _target_bed,
                '', // prefix for output files
                'g.vcf', // extension for the output files
                'HC_individually_genotyped_gvcf' // output directory name
                )

     ConcatVCF.out.concatenated_vcf_with_index
        .map { vc, idPatient, idSample, gvcfFile, tbiFile ->
            status = status_map[idPatient, idSample]
            gender = gender_map[idPatient]
            gvcf = "${params.outdir}/VariantCalling/${idSample}/HC_individually_genotyped_gvcf/${idSample}.g.vcf.gz"
            tbi = "${params.outdir}/VariantCalling/${idSample}/HC_individually_genotyped_gvcf/${idSample}.g.vcf.gz.tbi"
            gvcf_file = file(gvcf)
            tbi_file = file(tbi)
            "${idPatient}\t${gender}\t${status}\t${idSample}\t${gvcf_file}\t${tbi_file}\n"
        }.collectFile(
            name: "individually_genotyped_gvcfs.tsv", sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
        )

    emit: 
        gvcf_HC = HaplotypeCaller.out.gvcf_HC
        gvcf_GenotypeGVCFs = HaplotypeCaller.out.gvcf_GenotypeGVCFs
        gvcf_per_sample = ConcatVCF.out.concatenated_vcf_with_index
} //  end of wf_haplotypecaller







// STEP GATK HAPLOTYPECALLER.1

process HaplotypeCaller {
    label 'container_llab'
    label 'memory_singleCPU_task_sq'
    label 'cpus_8'
    
    tag {idSample + "-" + intervalBed.baseName}
    // tag {idSample} 
    // publishDir "${params.outdir}/VariantCalling/${idSample}/HaplotypeCaller", mode: params.publish_dir_mode
    input:
        tuple idPatient, idSample, file(bam), file(bai), file(intervalBed) 
        // tuple idPatient, idSample, file(bam), file(bai) 
        file(fasta)
        file(fastaFai)
        file(dict)
        file(dbsnp)
        file(dbsnpIndex)

    output:
        tuple val("HaplotypeCallerGVCF"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.g.vcf"), emit: gvcf_HC
        // tuple idPatient, idSample, file("${intervalBed.baseName}_${idSample}.g.vcf"), emit: gvcf_HaplotypeCaller
        tuple idPatient, idSample, file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf"), emit: gvcf_GenotypeGVCFs
        // tuple val("${intervalBed.baseName}"), idPatient, idSample, file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf"), emit: gvcf_GenotypeGVCFs
        

    when: 'haplotypecaller' in tools

    script:
    """
    init.sh
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -L ${intervalBed} \
        -D ${dbsnp} \
        -O ${intervalBed.baseName}_${idSample}.g.vcf \
        -ERC GVCF
    """
}

// process GvcfToVcf{
//     label 'memory_singleCPU_task_sq'
//     label 'cpus_2'
//     // label 'memory_max'
//     // label 'cpus_max'

//     tag {idSample + "-" + gvcf.baseName}
//     // tag {idSample} 
//     publishDir "${params.outdir}/VariantCalling/${idSample}/HC_individually_genotyped_vcf", mode: params.publish_dir_mode
//     input:
//         tuple val(variantCaller), idPatient, idSample, file(gvcf)

//     output:
//         // tuple val('HaplotypeCaller_Individually_Genotyped'), idPatient, idSample, file("${gvcf.simpleName}.vcf"), emit: vcf_HaplotypeCaller
//         tuple val('HaplotypeCaller_Individually_Genotyped'), idPatient, idSample, file(out_file), emit: vcf_HaplotypeCaller

//     when: 'haplotypecaller' in tools

//     script:
//     // fn=gvcf.fileName
//     // prefix=fn.minus(".g.vcf")
//     // out_file="${gvcf.fileName}.vcf"
//     prefix="${gvcf.fileName}" - ".g.vcf.gz"
//     // We'll first decompress the gvcf
//     in_file= "${gvcf.fileName}" - ".gz"
//     out_file="${prefix}.vcf"
//     """
//     gzip -d --force ${gvcf}
//     extract_variants < ${in_file} > ${out_file}
//     """
// }



// STEP GATK HAPLOTYPECALLER.1.5
