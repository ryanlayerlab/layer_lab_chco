tools = params.globals.tools

include {ConcatVCF} from './utility'
workflow wf_mutect2_single{
    take: _int_bam_recal
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _germline_resource
    take: _germline_resource_index
    take: _target_bed
    main:
        Mutect2Single(
            _int_bam_recal,
            _fasta,
            _fasta_fai,
            _dict,
            _germline_resource,
            _germline_resource_index
        )
        // group scattered vcf's (per interval) by the 
        // idPatient, and idSample

        _concat_vcf = Mutect2Single.out.vcf
                    .groupTuple(by: [0, 1])
                    
        _concat_vcf = Channel.from('Mutect2_single')
                    .combine(_concat_vcf)
                    //   .dump(tag: 'vcf_concat_vcf')
        
        ConcatVCF(
            _concat_vcf,
            _fasta_fai,
            _target_bed,
            'mutect2_single_unfiltered_', // prefix for output files
            'vcf', // extension for the output files
            'Mutect2_single_mode_unfiltered' // output directory name
            )
        // Now merge the stats together (groupby [idPatient, idSample])
        _merge_stats = Mutect2Single.out.stats
                    .groupTuple(by: [0, 1])

        MergeMutect2SingleStats(_merge_stats)
        // Now group per sample vcfs with correcsponding stats
        // remove the variant caller lable
        _concatenated_vcf = ConcatVCF.out.concatenated_vcf_with_index
                            .map{variantCaller, idPatient, idSample, vcf, tbi ->
                            [idPatient, idSample, vcf, tbi]}
        // join operator will join vcfs with their corresponding stats on the 
        // matching key [idPatient, idSample]
        _vcf_for_filtering = _concatenated_vcf
                            .join(MergeMutect2SingleStats.out, by:[0,1])
        FilterMutect2SingleCalls(_vcf_for_filtering,
                                _fasta,
                                _fasta_fai,
                                _dict,
                                _germline_resource,
                                _germline_resource_index)
         FilterMutect2SingleCalls.out
        

    emit:
        vcf = FilterMutect2SingleCalls.out
} // end of wf_mutect2_single
process Mutect2Single{
    label 'container_llab'
    tag {idSample + "-" + intervalBed.baseName}
    label 'cpus_16'

    input:
        tuple idPatient, idSample, file(bam), file(bai), 
            file(intervalBed)
        file(fasta)
        file(fastaFai)
        file(dict)
        file(germlineResource)
        file(germlineResourceIndex)

    output:
        tuple idPatient, idSample, file(out_vcf), emit: vcf
        tuple idPatient, idSample, file(out_stats) , emit: stats
    
    when: 'mutect2_single' in tools
    
    script:
    out_vcf = "${intervalBed.baseName}_${idSample}.vcf"
    out_stats = "${intervalBed.baseName}_${idSample}.vcf.stats"
    """
    # max-mnp-distance is set to 0 to avoid a bug in 
    # next process GenomicsDbImport
    # See https://gatk.broadinstitute.org/hc/en-us/articles/360046224491-CreateSomaticPanelOfNormals-BETA-
    init.sh
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
      Mutect2 \
      -R ${fasta} \
      -I ${bam}  \
      -max-mnp-distance 0 \
      -L ${intervalBed} \
      --germline-resource ${germlineResource} \
      -O ${out_vcf}
    """
}

process MergeMutect2SingleStats {
    label 'container_llab'
    label 'cpus_16'
    tag {idSample}

    // publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Mutect2", mode: params.publishDirMode

    input:
        tuple idPatient, idSample, file(statsFiles)// the actual stats files

    output:
        tuple idPatient, idSample, file("${idSample}.vcf.gz.stats")

    when: 'mutect2_single' in tools

    script:     
      stats = statsFiles.collect{ "-stats ${it} " }.join(' ')
    """
    init.sh
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        MergeMutectStats \
        ${stats} \
        -O ${idSample}.vcf.gz.stats
    """
} // end of MergeMutect2SingleStats


process FilterMutect2SingleCalls {
    label 'container_llab'
    label 'cpus_1'

    tag {idSample}

    publishDir "${params.outdir}/VariantCalling/${idSample}/mutect2_single_filtered", mode: params.publish_dir_mode

    input:
        tuple   idPatient, 
                idSample, 
                file(unfiltered), file(unfilteredIndex),
                file("${idSample}.vcf.gz.stats")
        
        file(fasta)
        file(fastaFai)
        file(dict)
        file(germlineResource)
        file(germlineResourceIndex)
        // file(intervals) from ch_intervals
        
    output:
        tuple val("Mutect2Single"), idPatient, idSample,
            file("mutect2_single_filtered_${idSample}.vcf.gz"),
            file("mutect2_single_filtered_${idSample}.vcf.gz.tbi"),
            file("mutect2_single_filtered_${idSample}.vcf.gz.filteringStats.tsv")

    // when: 'mutect2' in tools && params.pon
     when: 'mutect2_single' in tools

    script:
    """
    init.sh
    # do the actual filtering
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        FilterMutectCalls \
        -V ${unfiltered} \
        --stats ${idSample}.vcf.gz.stats \
        -R ${fasta} \
        -O mutect2_single_filtered_${idSample}.vcf.gz
    """
}