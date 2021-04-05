tools = params.globals.tools

workflow wf_gatk_cnv_somatic{
    take: _recal_bam
    take: _target_bed
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _read_count_somatic_pon
    
    main:
        /* GATK Somatic Copy Number related calls */
        /* Starting point is duplicated marked bams from MarkDuplicates.out.marked_bams with the following structure */
        /* MarkDuplicates.out.marked_bams => [idPatient, idSample, md.bam, md.bam.bai]*/

        PreprocessIntervals(_target_bed,
                            _fasta,
                            _fasta_fai,
                            _dict
                            )

        CollectReadCounts(_recal_bam,
                        PreprocessIntervals.out)

        // Create pon if specified
        // (normal_read_counts, tumor_read_counts) = 
        // CollectReadCounts.out
        // .branch{
        //      _: status_map[it[0], it[1]] == 0
        //     __: status_map[it[0], it[1]] == 1
        // }
        
        // normal_read_counts = 
        // CollectReadCounts.out.collect()
        // .filter{
        //      status_map[it[0], it[1]] == 0
        // }

        // normal_read_counts
        // .dump(tag: 'normal_read_counts: ')
        
        // normal_read_counts_hdf5 = 
        //         normal_read_counts
        //         .map{idPatient, idSample, hdf5 -> 
        //             [hdf5]
        //         }
        //         .dump(tag: 'hdf5: ')
        
        // tumor_read_counts
        // .dump(tag: 'tumor_read_counts: ')

        // CreateReadCountPon(
        //     normal_read_counts_hdf5.collect()
        // )
        
        DenoiseReadCounts(
            CollectReadCounts.out.sample_read_counts,
            _read_count_somatic_pon
        )

        PlotDenoisedCopyRatios(DenoiseReadCounts.out.denoised_cr,
                                _dict)

        denoised_cr_model_segments = DenoiseReadCounts.out.denoised_cr
                                    .map{ idPatien, idSample, std_cr, denoised_cr ->
                                        [idPatien, idSample, denoised_cr]
                                        }

        ModelSegments(denoised_cr_model_segments)
        plot_modeled_segments = ModelSegments.out.modeled_seg
                                .join(DenoiseReadCounts.out.denoised_cr, by: [0,1])
                                .map{idPatient, idSample, cr_seg, model_final_seg, std_cr, denoised_cr -> 
                                    [idPatient, idSample, model_final_seg, denoised_cr]
                                }
                                .dump(tag: 'plot_modeled_segments')

        PlotModeledSegments(plot_modeled_segments,
                            _dict)
        call_cr_seg = ModelSegments.out.modeled_seg                        
                                .map{idPatient, idSample, cr_seg, model_final_seg -> 
                                    [idPatient, idSample, cr_seg]
                                }
                                .dump(tag: 'call_cr_seg')

        CallCopyRatioSegments(call_cr_seg)
} // end of wf_gatk_somatic_cnv

/*
================================================================================
                                     Somatic CNV
================================================================================
*/

process PreprocessIntervals {
    label 'container_llab'
    label 'cpus_8'
    
    input:    
        file(intervalBed)
        file(fasta)
        file(fasta_fai)
        file(dict)
    
    output:
        // file("preprocessed_intervals.interval_list"), emit: 'processed_intervals'
        file("preprocessed_intervals.interval_list")

    // when: ('gatk_cnv_somatic' in tools) || ('gen_read_count_pon' in tools)
    when: ('gatk_cnv_somatic' in tools ||
           'gatk_cnv_germline_cohort_mode' in tools ||
           'gatk_cnv_germline' in tools
           )
    
    script:
    intervals_options = params.no_intervals ? "" : "-L ${intervalBed}"
    padding_options =  params.no_intervals ? "--padding 0" : "--padding 250"
    bin_options =  params.no_intervals ? "--bin-length 1000" : "--bin-length 0"

    """
    init.sh
    gatk PreprocessIntervals \
        ${intervals_options} \
        ${padding_options} \
        ${bin_options} \
        -R ${fasta} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O preprocessed_intervals.interval_list
    """
}

process CollectReadCounts {
    label 'container_llab'
    label 'cpus_32'
    tag "${idSample}"
    
    input:
        tuple idPatient, idSample, file(bam), file(bai)
        file(preprocessed_intervals)

    output:
        tuple idPatient, idSample, file("${idSample}.counts.hdf5"), emit: 'sample_read_counts'

    when: ('gatk_cnv_somatic' in tools ||
           'gatk_cnv_germline_cohort_mode' in tools ||
           'gatk_cnv_germline' in tools
           )

    script:
    """
    init.sh
    gatk CollectReadCounts \
        -I ${bam} \
        -L ${preprocessed_intervals} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ${idSample}.counts.hdf5
    """
}

process CreateReadCountPon {
    label 'container_llab'
    // echo true
    tag "ReadCountPon"
    
    publishDir "${params.outdir}/Preprocessing/ReadCountPon/", 
    mode: params.publish_dir_mode

    
    input:
    file(read_count_hdf5s)
    
    // file(this_read)

    output:
    file(out_file)

    script:
    when:'gen_read_count_pon' in tools
    // sample = this_read.simpleName
    out_file = "read_count_pon.hdf5"
    params_str = ''
    // Only get the normal samples
    read_count_hdf5s.each{
        params_str = "${params_str} -I ${it}"
    }

    
    """
    init.sh
    gatk CreateReadCountPanelOfNormals \
        $params_str \
        -O $out_file
    """
}

process DenoiseReadCounts {
    label 'container_llab'
    label 'cpus_32'
    tag "${idSample}"
    
    publishDir "${params.outdir}/Preprocessing/${idSample}/DenoisedReadCounts/", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file( "${idSample}.counts.hdf5")
        file(read_count_somatic_pon)

    output:
        tuple idPatient, idSample, file(std_copy_ratio), file(denoised_copy_ratio), emit: 'denoised_cr'

    when: ('gatk_cnv_somatic' in tools)

    script:
    std_copy_ratio = "${idSample}.standardizedCR.tsv"
    denoised_copy_ratio = "${idSample}.denoisedCR.tsv"
    pon_option = params.read_count_pon ? "--count-panel-of-normals ${read_count_somatic_pon}" : ""
    """
    init.sh
    gatk DenoiseReadCounts \
        -I ${idSample}.counts.hdf5 \
        ${pon_option} \
        --standardized-copy-ratios ${std_copy_ratio} \
        --denoised-copy-ratios ${denoised_copy_ratio}
    """
}

process PlotDenoisedCopyRatios {
    label 'container_llab'
    label 'cpus_16'
    tag "${idSample}"
    
    publishDir "${params.outdir}/Preprocessing/${idSample}/", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file(std_copy_ratio), file(denoised_copy_ratio)
        file(dict)
    
    output:
        file(out_dir)  

    when: ('gatk_cnv_somatic' in tools)

    script:
    out_dir = "PlotDenoisedReadCounts" 

    """
    init.sh
    mkdir ${out_dir}
    gatk PlotDenoisedCopyRatios \
        --standardized-copy-ratios ${std_copy_ratio} \
        --denoised-copy-ratios ${denoised_copy_ratio} \
        --sequence-dictionary ${dict} \
        --output-prefix ${idSample} \
        -O ${out_dir}
    """
}

process ModelSegments {
    label 'container_llab'
    label 'cpus_32'
    tag "${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}", mode: params.publish_dir_mode
    
    input:
         tuple idPatient, idSample, file(denoised_copy_ratio)

    output:
        tuple idPatient, idSample, file("${out_dir}/${idSample}.cr.seg"), file("${out_dir}/${idSample}.modelFinal.seg"), emit: 'modeled_seg'

    when: ('gatk_cnv_somatic' in tools)

    script:
    out_dir = "ModeledSegments"

    """
    init.sh
    mkdir $out_dir
    gatk ModelSegments \
        --denoised-copy-ratios ${denoised_copy_ratio} \
        --output-prefix ${idSample} \
        -O ${out_dir}
    """
}

process PlotModeledSegments {
    label 'container_llab'
    label 'cpus_8'
    tag "${idSample}"
    
    publishDir "${params.outdir}/VariantCalling/${idSample}", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file("${idSample}.modelFinal.seg"), file("${idSample}.denoisedCR.tsv")
        file(dict)
    output:
    file(out_dir)
    
    when: ('gatk_cnv_somatic' in tools)
    script:
    out_dir = "PlotsModeledSegments"
    
    """
    init.sh
    mkdir $out_dir
    gatk PlotModeledSegments \
        --denoised-copy-ratios ${idSample}.denoisedCR.tsv \
        --segments ${idSample}.modelFinal.seg \
        --sequence-dictionary ${dict} \
        --output-prefix ${idSample} \
        -O $out_dir
    """
}

process CallCopyRatioSegments {
    label 'container_llab'
   label 'cpus_8'
    tag "${idSample}"
    
    publishDir "${params.outdir}/VariantCalling/${idSample}/CalledCopyRatioSegments", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file("${idSample}.cr.seg")
    
    output:
        file("${idSample}.called.seg")

    when: ('gatk_cnv_somatic' in tools)
    script:
    
    """
    init.sh
    gatk CallCopyRatioSegments \
        -I ${idSample}.cr.seg \
        -O ${idSample}.called.seg
    """
}