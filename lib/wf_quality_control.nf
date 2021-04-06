process exonCoverage{
    label 'cpus_16'
    tag {idPatient + "-" + idSample}
    label 'container_llab'

    publishDir "${params.outdir}/Reports/${idSample}/exonCoverage/", mode: params.publish_dir_mode
    publishDir "${params.outdir}/QC/${idSample}/exonCoverage", mode: params.publish_dir_mode

    cache false

    input:
    tuple idPatient, idSample, file(bam), file(bai)
    // tuple idPatient, idSample, file(bam)
    file(fasta)
    file(fastaFai)
    file(dict)
    file(targetBED)
    file(baitBED)
    val outname


    output:
    path "${idSample}.*", emit: files
    file("target.interval_list")

    when: ! ('chco_qc' in _skip_qc)  && params.bait_bed

    script:
    """
    init.sh
    gatk BedToIntervalList -I ${targetBED} -O target.interval_list -SD ${dict}
    gatk BedToIntervalList -I ${baitBED} -O bait.interval_list -SD ${dict}
    echo "Hi"
    gatk --java-options -Xmx32G CollectHsMetrics --VALIDATION_STRINGENCY SILENT \
    -I ${bam} \
    -O ${idSample}.${outname}.hs_metrics.txt \
    -TI target.interval_list \
    -BI bait.interval_list \
    --PER_BASE_COVERAGE ${bam.baseName}.per_base_coverage.txt \
    -R ${fasta}
    exonCoverage.py target.interval_list ${bam.baseName}.per_base_coverage.txt ${bam.baseName}_per_exon_coverage.txt
    """
}

process onTarget{
    label 'container_llab'
    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Reports/${idSample}/onTarget/", mode: params.publish_dir_mode
    publishDir "${params.outdir}/QC/${idSample}/onTarget/", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(bam), file(bai)
    // tuple idPatient, idSample, file(bam)
    file(fasta)
    file(fastaFai)
    file(dict)
    file(probes)
    file(probes250)


    output:
    file("${bam.baseName}_on_target.txt")

    when: ! ('chco_qc' in _skip_qc)

    script:
    """
    init.sh
    gatk BedToIntervalList -I ${probes} -O probes.interval_list -SD ${dict}
    gatk BedToIntervalList -I ${probes250} -O probes250.interval_list -SD ${dict}
        bedtools intersect -a $bam -b ${probes250} | gatk --java-options -Xmx32G CollectHsMetrics \
    --VALIDATION_STRINGENCY SILENT \
        -I /dev/stdin \
        -O ${bam.baseName}_on_target.txt \
        -TI probes250.interval_list \
        -BI probes.interval_list \
        -R $fasta
    """
}

workflow wf_raw_bam_exonCoverage{

    take: _bams // tuple idPatient, idSample, file(bam), file(bai)
    take: _fasta
    take: _fasta_fai
    take: _dict
    take: _target
    take: _bait


    main:
        exonCoverage(_bams,_fasta,_fasta_fai,_dict,_target,_bait,"raw")

    emit:
        raw_onTarget = exonCoverage.out.files
}

workflow wf_qc_fingerprinting_sites{

    take: _bam
    take: _sites

    main:
         dnaFingerprint(_bam,_sites,"Extra")

    emit:
        fingerprint = dnaFingerprint.out
}

process insertSize{
    label 'container_llab'
    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Reports/${idSample}/insertSize/", mode: params.publish_dir_mode
    publishDir "${params.outdir}/QC/${idSample}/insertSize", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(bam), file(bai)


    output:
    path "${idSample}_insert_size_metrics.txt", emit: files
    file("${idSample}_insert_size_histogram.pdf")

    when: ! ('chco_qc' in _skip_qc)

    script:
    """
        gatk --java-options -Xmx32G CollectInsertSizeMetrics \
        -I $bam \
        -O ${idSample}_insert_size_metrics.txt \
        -H ${idSample}_insert_size_histogram.pdf \
        -M 0.5
    """
}

process dnaFingerprint{
    tag {idPatient + "-" + idSample}
    label 'container_llab'
    publishDir "${params.outdir}/Reports/${idSample}/FingerPrinting/${type}/", mode: params.publish_dir_mode
    publishDir "${params.outdir}/QC/${idSample}/FingerPrinting/${type}/", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(bam), file(bai)
    file(finger_printing_sites)
    val type


    output:
    file("${idSample}_DNA_Fingerprint.txt")

    when: ! ('chco_qc' in _skip_qc)

    script:
    """
        dnaFingerPrinting.py $bam $finger_printing_sites $idSample
    """
}

process collectQC{
    label 'container_py3_pandas'
    publishDir "${params.outdir}/Reports/${idSample}/FingerPrinting/", mode: params.publish_dir_mode
    publishDir "${params.outdir}/QC/collectQC", mode: params.publish_dir_mode

    cache false

    input:
    file(sample_file)
    file(results_dir)
    file(exon)
    file(raw_exon)
    file(insertsize)
    file(fingerprint)
    file(bcf)
    file(unknown)

    when: ! ('chco_qc' in _skip_qc)

    output:
    file("QC_Stats.xlsx")


    script:
    """
        echo ${sample_file}
        echo ${params.outdir}
        echo ${params.input}
        echo \$PWD
        collectQC.py ${sample_file} ${params.outdir} \$PWD
    """
}

process add_somalier_to_QC{
    label 'container_py3_pandas'
    publishDir "${params.outdir}/QC/collectQC", mode: params.publish_dir_mode

    input:
    tuple file(html), file(pairs), file(samples)
    file(pedigree)
	file(pre_QC_stats)


    when: ! ('chco_qc' in _skip_qc)

    output:
    file("QC_Stats_Final.xlsx")


    script:
    """
        somalier_to_excel.py $pre_QC_stats $samples $pairs $pedigree
    """
}

process add_cohort_vc_to_qc_report{
    tag {idPatient + "-" + idSample}
    label 'container_py3_pandas'

    publishDir "${params.outdir}/QC/collectQC", mode: params.publish_dir_mode

    input:
    tuple file(vcfgz), file(vcfgzindex)
    file(qc_file)

    output:
    file('QC_Stats_*.xlsx')

    script:
    """
    zcat $vcfgz | add_sample_count_to_cohort_vcf.py > cohort_vcf_with_count_column.tsv
    # add it to the QC report now
    add_cohort_vcf_to_qc.py $qc_file cohort_vcf_with_count_column.tsv Cohort_VCF
    """
}

process add_cohort_CNVs_to_qc_report{
    tag {idPatient + "-" + idSample}
    label 'container_py3_pandas'

    publishDir "${params.outdir}/QC/collectQC", mode: params.publish_dir_mode

    input:
    file(vcf)
    file(qc_file)

    output:
    file('QC_Stats*.xlsx')

    script:
    """
    cat $vcf | add_sample_count_to_cohort_vcf.py > cohort_vcf_with_count_column.tsv
    # add it to the QC report now
    add_cohort_vcf_to_qc.py $qc_file cohort_vcf_with_count_column.tsv Merged_CNV
    """
}
