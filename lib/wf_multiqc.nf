/*
================================================================================
                                     MultiQC
================================================================================
*/

// STEP MULTIQC

skip_qc = params.globals.skip_qc

workflow wf_multiqc{
    // to do
    take: software_versions
    take: bam_qc
    take: fastqc_out 
    take: bcftools_out 
    take: vcftools_out 
    // take: dm_bam_stats
    take: samtools_stats
    take: alignment_summary_metrics
    take: insert_size_metrics
    take: hs_metrics

    main:
        if (!('multiqc' in skip_qc)){
            MultiQC(
            Channel.value(""),
            software_versions,
            bam_qc,
            fastqc_out,
            bcftools_out,
            vcftools_out,
            // dm_bam_stats,
            samtools_stats,
            alignment_summary_metrics,
            insert_size_metrics,
            hs_metrics
            // SnpEff.out.snpEff_report,
            )
        }
        
} // end of wf_germline_cnv


process MultiQC {
    label 'container_sarek'
    memory '64 GB'
    cpus 8    

    publishDir "${params.outdir}/Reports/MultiQC", mode: params.publish_dir_mode
    input:
        file (multiqcConfig) 
        file (versions) 
        file ('bamQC/*') 
        file ('FastQC/*') 
        file ('BCFToolsStats/*') 
        file ('VCFTools/*')
        // file ('MarkDuplicates/*') 
        file ('SamToolsStats/*') 
        file ('CollectAlignmentSummary/*')
        file ('CollectInsertSizeMetrics/*')
        file ('CollectHsMetrics/*')
        // file ('snpEff/*') 

    output:
        tuple file("*multiqc_report.html"), file("*multiqc_data") 

    script:
    """
    multiqc -f -v .
    """
}
