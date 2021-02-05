tools = params.globals.tools
include {ConcatVCF} from './utility'
workflow wf_mutect2_TN{
    take: _int_pair_bam
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _germline_resource
    take: _germline_resource_index
    take: _pon_somatic
    take: _pon_somatic_index
    take: _target_bed

    main:
        Mutect2TN(
            _int_pair_bam,
            _fasta,
            _fasta_fai,
            _dict,
            _germline_resource,
            _germline_resource_index,
            _pon_somatic,
            _pon_somatic_index
        )
        // group scattered vcf's (per interval) by the 
        // idPatient, and idSample

        _concat_vcf = Mutect2TN.out.vcf
                      .groupTuple(by: [0, 1])
                      
        _concat_vcf = Channel.from('Mutect2_TN')
                      .combine(_concat_vcf)
                    //   .dump(tag: 'vcf_concat_vcf')
        
        ConcatVCF(
            _concat_vcf,
            _fasta_fai,
            _target_bed,
            'mutect2_TN_unfiltered_', // prefix for output files
            'vcf', // extension for the output files
            'Mutect2_unfiltered' // output directory name
            )
        // Now merge the stats together (groupby [idPatient, idSample])
         _merge_stats = Mutect2TN.out.stats
                      .groupTuple(by: [0, 1, 2])

        MergeMutect2TNStats(_merge_stats)
        // Now group per sample vcfs with correcsponding stats
        // remove the variant caller lable
        _concatenated_vcf = ConcatVCF.out.concatenated_vcf_with_index
                            .map{variantCaller, idPatient, idSampleTN, vcf, tbi ->
                            [idPatient, idSampleTN, vcf, tbi]}
        // join operator will join vcfs with their corresponding stats on the 
        // matching key [idPatient, idSample]
        _vcf_for_filtering = _concatenated_vcf
                            .join(MergeMutect2TNStats.out, by:[0,1])
        FilterMutect2TNCalls(_vcf_for_filtering,
                                 _fasta,
                                 _fasta_fai,
                                 _dict,
                                 _germline_resource,
                                 _germline_resource_index)

    emit:
        vcf = FilterMutect2TNCalls.out
} // end of wf_mutect2_TN
/*
================================================================================
                                     Mutect2 (Tumor/Normal Mode)
================================================================================
*/

process Mutect2TN{
    label 'container_llab'
    tag {idSampleTumor + "_vs_" + idSampleNormal + "-" + intervalBed.baseName}
    label 'cpus_2'

    input:
        tuple idPatient, 
            idSampleNormal, file(bamNormal), file(baiNormal),
            idSampleTumor, file(bamTumor), file(baiTumor), 
            file(intervalBed)
        file(fasta)
        file(fastaFai)
        file(dict)
        file(germlineResource)
        file(germlineResourceIndex)
        file(ponSomatic)
        file(ponSomaticIndex)

    output:
        tuple idPatient,
            val("${idSampleTumor}_vs_${idSampleNormal}"),
            file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf"), emit: vcf
        
        tuple idPatient,
            idSampleTumor,
            idSampleNormal,
            file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf.stats"), emit: stats

    when: 'mutect2' in tools

    script:
    // please make a panel-of-normals, using at least 40 samples
    // https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
    PON = params.somatic_pon ? "--panel-of-normals ${ponSomatic}" : ""
    // PON =  ""
    """
    init.sh
    # Get raw calls
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
      Mutect2 \
      -R ${fasta}\
      -I ${bamTumor}  -tumor ${idSampleTumor} \
      -I ${bamNormal} -normal ${idSampleNormal} \
      -L ${intervalBed} \
      --germline-resource ${germlineResource} \
      ${PON} \
      -O ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
    """
}

// // STEP GATK MUTECT2.2 - MERGING STATS

process MergeMutect2TNStats {
    label 'container_llab'
    label 'cpus_16'
    tag {idSampleTumor + "_vs_" + idSampleNormal}

    publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Mutect2", mode: params.publish_dir_mode

    input:
        // tuple caller, idPatient, idSampleTumor_vs_idSampleNormal, file(vcfFiles) // corresponding small VCF chunks
        tuple idPatient, idSampleTumor, idSampleNormal, file(statsFiles)// the actual stats files

    output:
        tuple idPatient,
            val("${idSampleTumor}_vs_${idSampleNormal}"),
            file("${idSampleTumor}_vs_${idSampleNormal}.vcf.gz.stats")

    when: 'mutect2' in tools

    script:     
      stats = statsFiles.collect{ "-stats ${it} " }.join(' ')
    """
    init.sh
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        MergeMutectStats \
        ${stats} \
        -O ${idSampleTumor}_vs_${idSampleNormal}.vcf.gz.stats
    """
} // end of MergeMutect2Stats

// // STEP GATK MUTECT2.3 - GENERATING PILEUP SUMMARIES

// process PileupSummariesForMutect2 {
//     tag {idSampleTumor + "_vs_" + idSampleNormal + "_" + intervalBed.baseName }
//     label 'cpus_1'

//     input:
//         tuple idPatient, 
//             idSampleNormal, file(bamNormal), file(baiNormal), 
//             idSampleTumor, file(bamTumor), file(baiTumor), 
//             file(intervalBed)
        
//         // set idPatient, idSampleNormal, idSampleTumor, file(statsFile) from intervalStatsFiles
//         file(germlineResource)
//         file(germlineResourceIndex)

//     output:
//         tuple idPatient,
//             idSampleTumor,
//             file("${intervalBed.baseName}_${idSampleTumor}_pileupsummaries.table")

//     when: 'mutect2' in tools && params.pon

//     script:
//     """
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         GetPileupSummaries \
//         -I ${bamTumor} \
//         -V ${germlineResource} \
//         -L ${intervalBed} \
//         -O ${intervalBed.baseName}_${idSampleTumor}_pileupsummaries.table
//     """
// } // end of PileupSummariesForMutect2

// // STEP GATK MUTECT2.4 - MERGING PILEUP SUMMARIES

// process MergePileupSummaries {
//     label 'cpus_1'

//     tag {idPatient + "_" + idSampleTumor}

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}/Mutect2", mode: params.publishDirMode

//     input:
//         tuple idPatient, idSampleTumor, file(pileupSums)
//         file(dict)

//     output:
//         file("${idSampleTumor}_pileupsummaries.table.tsv")

//     when: 'mutect2' in tools
//     script:
//         allPileups = pileupSums.collect{ "-I ${it} " }.join(' ')
//     """
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         GatherPileupSummaries \
//         --sequence-dictionary ${dict} \
//         ${allPileups} \
//         -O ${idSampleTumor}_pileupsummaries.table.tsv
//     """
// }

// // STEP GATK MUTECT2.5 - CALCULATING CONTAMINATION

// process CalculateContamination {
//     label 'cpus_1'

//     tag {idSampleTumor + "_vs_" + idSampleNormal}

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}/Mutect2", mode: params.publishDirMode

//     input:
//         // tuple idPatient, 
//         //     idSampleNormal, file(bamNormal), file(baiNormal), 
//         //     idSampleTumor, file(bamTumor), file(baiTumor) from pairBamCalculateContamination 
//         file("${idSampleTumor}_pileupsummaries.table")
  
//     output:
//         file("${idSampleTumor}_contamination.table") into contaminationTable

//     when: 'mutect2' in tools && params.pon

//     script:     
//     """
//     # calculate contamination
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         CalculateContamination \
//         -I ${idSampleTumor}_pileupsummaries.table \
//         -O ${idSampleTumor}_contamination.table
//     """
// }

// // STEP GATK MUTECT2.6 - FILTERING CALLS

process FilterMutect2TNCalls {
    label 'container_llab'
    label 'cpus_1'

    tag {idSampleTN}

    publishDir "${params.outdir}/VariantCalling/${idSampleTN}/Mutect2", mode: params.publish_dir_mode

    input:
        // tuple variantCaller, 
        //     idPatient, idSampleTN, file(unfiltered), file(unfilteredIndex),
        //     file("${idSampleTN}.vcf.gz.stats"),
        //     file("${idSampleTN}_contamination.table")
        tuple idPatient, 
            idSampleTN, 
            file(unfiltered), file(unfilteredIndex),
            file("${idSampleTN}.vcf.gz.stats")
            // file("${idSampleTN}_contamination.table")
        
        file(fasta)
        file(fastaFai)
        file(dict)
        file(germlineResource)
        file(germlineResourceIndex)
        // file(intervals) from ch_intervals
        
    output:
        tuple val("Mutect2"), idPatient, idSampleTN,
            file("filtered_mutect2_${idSampleTN}.vcf.gz"),
            file("filtered_mutect2_${idSampleTN}.vcf.gz.tbi"),
            file("filtered_mutect2_${idSampleTN}.vcf.gz.filteringStats.tsv")

    // when: 'mutect2' in tools && params.pon
    when: 'mutect2' in tools

    script:
    """
    init.sh
    # do the actual filtering
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        FilterMutectCalls \
        -V ${unfiltered} \
        --stats ${idSampleTN}.vcf.gz.stats \
        -R ${fasta} \
        -O filtered_mutect2_${idSampleTN}.vcf.gz
    """
}