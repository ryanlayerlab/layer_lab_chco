tools = params.globals.tools
include { PrepareTargets; 
          PrepareAntiTargets; 
          SampleAntiTargetCoverage;
          SampleTargetCoverage
        } from './wf_cnvkit_single'

workflow wf_cnvkit_gen_ref{
    take: _md_bam_collected
    take: _cnv_target_bed
    take: _fasta
    take: _fasta_fai
    main:
     /* CNVKit Somatic Copy Number related calls */
    /* Starting point is duplicated marked bams from MarkDuplicates.out.marked_bams with the following structure */
    /* MarkDuplicates.out.marked_bams => [idPatient, idSample, md.bam, md.bam.bai]*/
        // CNVKitSingle(
        //     _md_bam_collected,
        //     _fasta,
        //     _target_bed)
        PrepareTargets(_cnv_target_bed)
        PrepareAntiTargets(_cnv_target_bed)
        _cnvkit_targets = PrepareTargets.out
        _cnvkit_antitargets = PrepareAntiTargets.out
        SampleTargetCoverage(_md_bam_collected,
                             _cnvkit_targets)

        SampleAntiTargetCoverage(_md_bam_collected,
                             _cnvkit_antitargets)
        _all_targ_cov = SampleTargetCoverage.out.collect()
        // _all_targ_cov.dump(tag: "all_target_coverages: ")
        
        _all_antitarg_cov = SampleAntiTargetCoverage.out.collect()
        ch_all_coverages = _all_targ_cov.mix(_all_antitarg_cov)
                            
        GenRef(_fasta,
                _fasta_fai,
                _all_targ_cov,
                _all_antitarg_cov
                )

} // end of wf_germline_cnv

process GenRef{
    label 'cpus_16'
    label 'container_llab'
    publishDir "${params.outdir}/Preprocessing/CNVKitRef/", mode: params.publish_dir_mode
    
    input:
        file(fasta)
        file(fasta_fai)
        file("tar_covs/*")
        file("anti_covs/*")


    output:
    file("Reference.cnn")

    when: 'cnvkit_gen_ref' in tools

    script:
    
    """
    init.sh
    cnvkit.py reference {tar_covs,anti_covs}/*coverage.cnn -f ${fasta} -o Reference.cnn
    """
}