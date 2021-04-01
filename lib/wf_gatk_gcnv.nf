tools =  params.globals.tools


// More infor at
// https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variant

// include {PreprocessIntervals} from './wf_gatk_cnv_somatic' 
// include {CollectReadCounts} from './wf_gatk_cnv_somatic' 
workflow wf_gatk_gcnv{
    take: _in_bam
    take: _target_bed
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _contig_ploidy_priors
    take: _ploidy_model
    take: _cohort_model
    
    main:
        /* GATK Somatic Copy Number related calls */
        /* Starting point is duplicated marked bams from MarkDuplicates.out.marked_bams with the following structure */
        /* MarkDuplicates.out.marked_bams => [idPatient, idSample, md.bam, md.bam.bai]*/
        // _in_bam.dump(tag:'_in_bam: ')
        // _target_bed.view{"_target_bed: $it"}
        // _contig_ploidy_priors.dump(tag:'_contig_ploidy_priors: ')
        // _ploidy_model.dump(tag:'_ploidy_model: ')
        // _cohort_model.dump(tag:'_cohort_model: ')
        // TestProc()
        PreprocessIntervals(_target_bed,
                            _fasta,
                            _fasta_fai,
                            _dict
                            )
        preprocessed_intervals = PreprocessIntervals.out
        CollectReadCounts(_in_bam,
                        preprocessed_intervals)

        AnnotateIntervals(_fasta,
                          _fasta_fai,
                          _dict,
                          preprocessed_intervals
                          )
        anno_intervals = AnnotateIntervals.out
        all_cvgs = CollectReadCounts.out.collect()
                    // .dump(tag: 'all_cvgs: ')
        FilterIntervals(preprocessed_intervals,
                        anno_intervals,
                        all_cvgs
                        )
        filtered_intervals = FilterIntervals.out
        // ploidy_calls = Channel.empty()
        // ploidy_model = Channel.empty()
        // g_cnv_model = Channel.empty()
        if('gatk_gcnv_cohort_mode' in tools){
            if (! params.gcnv_contig_ploidy_priors)
                error "Running 'gatk_gcnv__cohort_mode' requires 'gcnv_contig_ploidy_priors'"
            
            DetermineGermlineContigPloidyCohortMode(
                filtered_intervals,
                _contig_ploidy_priors,
                all_cvgs
            )
            ploidy_calls = DetermineGermlineContigPloidyCohortMode.out[0]
            //ploidy_model = DetermineGermlineContigPloidyCohortMode.out[1]
            GermlineCNVCallerCohortMode(
                anno_intervals,
                filtered_intervals,
                ploidy_calls,
                all_cvgs
            )
            //g_cnv_model = GermlineCNVCallerCohortMode()

        }
        
        if('gatk_gcnv' in tools){
            if (! params.gcnv_ploidy_model)
                error "Running 'gatk_gcnv' requires 'gcnv_ploidy_model'"
            if (! params.gcnv_cohort_model)
                error "Running 'gatk_gcnv' requires 'gcnv_cohort_model'"
        }

        DetermineGermlineContigPloidyCaseMode(
            CollectReadCounts.out,
            _ploidy_model
        )
        cvg_plus_ploidy = CollectReadCounts.out
                        .join(DetermineGermlineContigPloidyCaseMode.out, by:[0,1])
                        .dump(tag: 'cvg_plus_ploidy: ')
        
        GermlineCNVCallerCaseMode(
            cvg_plus_ploidy,
            _cohort_model,
            filtered_intervals
        )
} // end of wf_gatk_somatic_cnv

/*
================================================================================
                                     germline CNV
================================================================================
*/
// process TestProc{
//     label 'container_llab'
//     label 'cpus_2'
//     echo true
//     // when: ('gatk_gcnv_cohort_mode' in tools )
//     // when: ('gatk_gcnv_cohort_mode' in tools )
//     script:
//     """
//     echo "TestProc says tools: ${tools}"
//     """

// }
process PreprocessIntervals {
    label 'container_llab'
    label 'cpus_8'
    publishDir "${params.outdir}/Preprocessing/gatk_gcnv/", mode: params.publish_dir_mode

    input:    
        file(intervalBed)
        file(fasta)
        file(fasta_fai)
        file(dict)
    
    output:
        // file("preprocessed_intervals.interval_list"), emit: 'processed_intervals'
        file("preprocessed_intervals.interval_list")

    // when: ('gatk_cnv_somatic' in tools) || ('gen_read_count_pon' in tools)
    when: ('gatk_gcnv_cohort_mode' in tools ) ||
           ('gatk_gcnv' in tools )
    
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
    publishDir "${params.outdir}/Preprocessing/gatk_gcnv/read_counts", mode: params.publish_dir_mode
    input:
        tuple idPatient, idSample, file(bam), file(bai)
        file(preprocessed_intervals)

    output:
        tuple idPatient, idSample, file("${idSample}.counts.tsv"), emit: 'sample_cvf'

    when: ('gatk_gcnv_cohort_mode' in tools ||
           'gatk_gcnv' in tools
           )

    script:
    """
    init.sh
    gatk CollectReadCounts \
        -I ${bam} \
        -L ${preprocessed_intervals} \
        --interval-merging-rule OVERLAPPING_ONLY \
        --format TSV \
        -O ${idSample}.counts.tsv
    """
}

process AnnotateIntervals {
    label 'container_llab'
    label 'cpus_16'
    publishDir "${params.outdir}/Preprocessing/gatk_gcnv/", mode: params.publish_dir_mode
    input:
        file(fasta)
        file(fasta_fai)
        file(dict)
        file(preprocessed_intervals)

    output:
        file("preprocessed_intervals.annotated.tsv")

    when: ('gatk_gcnv_cohort_mode' in tools)

    script:
    """
    init.sh
     gatk AnnotateIntervals \
        -L ${preprocessed_intervals} \
        -R ${fasta} \
        -imr OVERLAPPING_ONLY \
        -O preprocessed_intervals.annotated.tsv
    """
}

process FilterIntervals {
    label 'container_llab'
    label 'cpus_16'
    echo true
    publishDir "${params.outdir}/Preprocessing/gatk_gcnv/", mode: params.publish_dir_mode
    input:
        file(preprocessed_intervals)
        file(annotated_intervals)
        file("cvg/*")

    output:
        file("cohort.gc.filtered.interval_list")
        

    when: ('gatk_gcnv_cohort_mode' in tools)

    script:
    """
    init.sh
    #cvg_opts=''
    for x in `ls cvg/*.tsv`
    do
        echo -n "-I \$x " >> args.list
    done
    
    #echo "my cvg_opts: \$cvg_opts"
    
    gatk FilterIntervals \
        -L ${preprocessed_intervals} \
        --annotated-intervals ${annotated_intervals} \
        --arguments_file args.list \
        -imr OVERLAPPING_ONLY \
        -O cohort.gc.filtered.interval_list

    """
}
// Optional step. Just for cohort model, to generate the model

process DetermineGermlineContigPloidyCohortMode {
    // label 'container_llab'
    label 'container_gatk'
    label 'cpus_32'
    // echo true
    publishDir "${params.outdir}/Preprocessing/gatk_gcnv/", mode: params.publish_dir_mode
    input:
        file(filterd_intervals)
        file(all_contig_ploidy_priors)
        file("cvg/")
        

    output:
        file("ploidy-calls")
        file("ploidy-model")

    when:  'gatk_gcnv_cohort_mode' in tools 

    script:
    """
    init.sh
    for x in `ls cvg/*.tsv`
    do
        echo -n "-I \$x " >> args.list
    done

    gatk DetermineGermlineContigPloidy \
    -L ${filterd_intervals} \
    --interval-merging-rule OVERLAPPING_ONLY \
    --arguments_file args.list \
    --contig-ploidy-priors ${all_contig_ploidy_priors} \
    --output . \
    --output-prefix ploidy \
    --verbosity DEBUG

    """
}

process DetermineGermlineContigPloidyCaseMode {
    // label 'container_llab'
    label 'container_gatk'
    label 'cpus_32'
    
    publishDir "${params.outdir}/Preprocessing/${idSample}/GATK_gcnv/", mode: params.publish_dir_mode
    input:
        tuple idPatient, idSample, file(sample_cvg)
        file(ploidy_model)

    output:
        tuple idPatient, idSample, file("ploidy-case-calls")

    when: ( 'gatk_gcnv' in tools )

    script:
     
    """
    init.sh
    gatk DetermineGermlineContigPloidy \
    --model ${ploidy_model} \
        -I ${sample_cvg} \
        -O . \
        --output-prefix ploidy-case \
        --verbosity DEBUG
    """
}

process GermlineCNVCallerCohortMode {
    // label 'container_llab'
    label 'container_gatk'
    label 'cpus_32'
    publishDir "${params.outdir}/Preprocessing/GATK_gcnv/", mode: params.publish_dir_mode
    input:
        file(annotated_intervals)
        file(filtered_intervals)
        file(ploidy_calls)
        file("cvg/*")

    output:
        file("gcnv_model")

    when:  'gatk_gcnv_cohort_mode' in tools 

    script:
     
    """
    init.sh
    for x in `ls cvg/*.tsv`
    do
        echo -n "-I \$x " >> args.list
    done
    gatk GermlineCNVCaller \
        --run-mode COHORT \
        -L ${filtered_intervals} \
        --arguments_file args.list \
        --contig-ploidy-calls ${ploidy_calls} \
        --annotated-intervals ${annotated_intervals} \
        --interval-merging-rule OVERLAPPING_ONLY \
        --output gcnv_model \
        --output-prefix gcnv_model \
        --verbosity DEBUG
    """
}

process GermlineCNVCallerCaseMode {
    // label 'container_llab'
    label 'container_gatk'
    label 'cpus_32'
    publishDir "${params.outdir}/VariantCalling/${idSample}/gcnv/", mode: params.publish_dir_mode
    input:
        tuple idPatient, idSample, file(sample_cvg), file(case_ploidy_calls)
        file(gcnv_cohort_model)
        file(filtered_intervals)
        // file(ploidy_calls)
        

    output:
        file("gCNV_model")

    when: ( 'gatk_gcnv' in tools )

    script:
     
    """
    init.sh
     gatk GermlineCNVCaller \
        --run-mode CASE \
        -I ${sample_cvg} \
        --contig-ploidy-calls ${case_ploidy_calls} \
        --model ${gcnv_cohort_model} \
        --output "${idSample}_vs_cohort" \
        --output-prefix  ${idSample}_vs_cohort\
        --verbosity DEBUG
    """
}