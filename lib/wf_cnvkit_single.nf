tools = params.globals.tools
workflow wf_cnvkit_single{
    take: _md_bam_collected
    take: _target_bed
    take: _fasta
    take: _fasta_fai
    take: _cnvkit_ref
    main:
     /* CNVKit Somatic Copy Number related calls */
    /* Starting point is duplicated marked bams from MarkDuplicates.out.marked_bams with the following structure */
    /* MarkDuplicates.out.marked_bams => [idPatient, idSample, md.bam, md.bam.bai]*/
        // CNVKitSingle(
        //     _md_bam_collected,
        //     _fasta,
        //     _target_bed)
        PrepareTargets(_target_bed)
        PrepareAntiTargets(_target_bed)
        _cnvkit_targets = PrepareTargets.out
        _cnvkit_antitargets = PrepareAntiTargets.out
        _ref = _cnvkit_ref
        if (!params.cnvkit_ref){
            PrepareFlatReference(_fasta,
                                 _fasta_fai,
                                 _cnvkit_targets,
                                 _cnvkit_antitargets
                                 )
            _ref = PrepareFlatReference.out
        }
        SampleTargetCoverage(_md_bam_collected,
                             _cnvkit_targets)

        SampleAntiTargetCoverage(_md_bam_collected,
                             _cnvkit_antitargets)
        
        _ch_fix_for_biases = 
                        SampleTargetCoverage.out
                        .join(SampleAntiTargetCoverage.out, by: [0,1])
                        .dump(tag: "target_cov_joined_by_anti: ")
        FixForBiases(_ch_fix_for_biases,
                     _ref)
        Segment(FixForBiases.out[0])
        CallSegment(Segment.out[0])
        Scatter(Segment.out[1])
        // Diagram(Segment.out[1])

} // end of wf_germline_cnv

process PrepareTargets{
    label 'cpus_8'
    label 'container_llab'
    publishDir "${params.outdir}/VariantCalling/CNVKit_preprocessing", mode: params.publish_dir_mode
    
    input:
        // file(fasta)
        // file(fastaFai)
        file(targetBED)
    
    output:
    file("cnvkit_targets.bed")

    when: 'cnvkit_single' in tools

    script:
    
    """
    init.sh
    # Generate the targets.bed
    cnvkit.py target ${targetBED} --split -o cnvkit_targets_with_chr.bed
    # remove the chr from the chromosoms
    cat cnvkit_targets_with_chr.bed | sed 's/^chr//' > cnvkit_targets.bed
    """
}

process PrepareAntiTargets{
    label 'cpus_8'
    label 'container_llab'
    publishDir "${params.outdir}/VariantCalling/CNVKit_preprocessing", mode: params.publish_dir_mode
    
    input:
        // file(fasta)
        // file(fastaFai)
        file(targetBED)
    
    output:
    file("cnvkit_antitargets.bed")

    when: 'cnvkit_single' in tools

    script:
    
    """
    init.sh
    cnvkit.py antitarget ${targetBED}  -o cnvkit_antitargets_with_chr.bed
    cat cnvkit_antitargets_with_chr.bed | sed 's/^chr//' > cnvkit_antitargets.bed
    """
}

process PrepareFlatReference{
    label 'cpus_8'
    label 'container_llab'
    publishDir "${params.outdir}/VariantCalling/CNVKit_preprocessing", mode: params.publish_dir_mode
    
    input:
        file(fasta)
        file(fastaFai)
        file(cnvkit_target)
        file(cnvkit_antitarget)
    
    output:
    file("flat_reference.cnn")

    when: !params.cnvkit_ref && 'cnvkit_single' in tools

    script:
    
    """
    init.sh
    cnvkit.py reference -o flat_reference.cnn -f ${fasta} -t ${cnvkit_target} -a ${cnvkit_antitarget}
    """
}

process SampleTargetCoverage{
    label 'cpus_16'
    label 'container_llab'
    publishDir "${params.outdir}/VariantCalling/${idSample}/CNVKit", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file(bam), file(bai)
        // file(fasta)
        // file(fastaFai)
        file(cnvkit_targets)
    
    output:
    //tuple idPatient, idSample, file("${idSample}.targetcoverage.cnn"), file("${idSample}.antitargetcoverage.cnn")
    tuple idPatient, idSample, file("${idSample}.targetcoverage.cnn")

    when: 'cnvkit_single' in tools

    script:
    
    """
    init.sh
    cnvkit.py coverage ${bam} ${cnvkit_targets} -o ${idSample}.targetcoverage.cnn
    """
}

process SampleAntiTargetCoverage{
    label 'cpus_16'
    label 'container_llab'
    publishDir "${params.outdir}/VariantCalling/${idSample}/CNVKit", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file(bam), file(bai)
        // file(fasta)
        // file(fastaFai)
        file(cnvkit_antitargets)
    
    output:
    tuple idPatient, idSample, file("${idSample}.antitargetcoverage.cnn")

    when: 'cnvkit_single' in tools

    script:
    
    """
    init.sh
    cnvkit.py coverage ${bam} ${cnvkit_antitargets} -o ${idSample}.antitargetcoverage.cnn
    """
}

process FixForBiases{
    label 'cpus_16'
    label 'container_llab'
    publishDir "${params.outdir}/VariantCalling/${idSample}/CNVKit", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file(target_coverage), file(antitarget_coverage)
        file(cnvkit_ref)
    
    output:
    tuple idPatient, idSample, file("${idSample}.cnr")

    when: 'cnvkit_single' in tools

    script:
    
    """
    init.sh
    cnvkit.py fix ${target_coverage} ${antitarget_coverage} ${cnvkit_ref} -o ${idSample}.cnr
    """
}

process Segment{
    label 'cpus_16'
    label 'container_llab'
    publishDir "${params.outdir}/VariantCalling/${idSample}/CNVKit", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file(sample_cnr)
        
    
    output:
    tuple idPatient, idSample, file("${idSample}.cns")
    tuple idPatient, idSample, file(sample_cnr), file("${idSample}.cns")
    // tuple idPatient, idSample, file("${idSample}.cns"), emit: cns
    // tuple idPatient, idSample, file(sample_cnr), file("${idSample}.cns"), emit:cnr_and_cns

    when: 'cnvkit_single' in tools

    script:
    
    """
    init.sh
    cnvkit.py segment ${sample_cnr} -o ${idSample}.cns
    """
}

process CallSegment{
    label 'cpus_16'
    label 'container_llab'
    publishDir "${params.outdir}/VariantCalling/${idSample}/CNVKit", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file(sample_cns)
        
    
    output:
    tuple idPatient, idSample, file("${idSample}.called.cns")
    // tuple idPatient, idSample, file(sample_cnr), file("${idSample}.cns")
    // tuple idPatient, idSample, file("${idSample}.cns"), emit: cns
    // tuple idPatient, idSample, file(sample_cnr), file("${idSample}.cns"), emit:cnr_and_cns

    when: 'cnvkit_single' in tools

    script:
    
    """
    init.sh
    cnvkit.py call ${sample_cns} -o ${idSample}.called.cns
    """
}

process Scatter{
    label 'cpus_16'
    label 'container_llab'
    publishDir "${params.outdir}/VariantCalling/${idSample}/CNVKit", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file(sample_cnr), file(sample_cns) 
        
    
    output:
    tuple idPatient, idSample, file("${idSample}_scatter.pdf")

    when: 'cnvkit_single' in tools

    script:
    
    """
    init.sh
    cnvkit.py scatter ${sample_cnr} -s ${sample_cns} -o ${idSample}_scatter.pdf
    """
}

process Diagram{
    label 'cpus_16'
    label 'container_llab'
    publishDir "${params.outdir}/VariantCalling/${idSample}/CNVKit", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file(sample_cnr), file(sample_cns)
        
    
    output:
    tuple idPatient, idSample, file("${idSample}_diagram.pdf")

    when: 'cnvkit_single' in tools

    script:
    
    """
    init.sh
    cnvkit.py diagram ${sample_cnr} -s ${sample_cns} -o ${idSample}_diagram.pdf
    """
}

process CNVKitSingle{
    label 'cpus_8'
    label 'container_llab'
    publishDir "${params.outdir}/VariantCalling/${idSample}/CNVKit", mode: params.publish_dir_mode
    
    input:
        tuple idPatient, idSample, file(bam), file(bai)
        file(fasta)
        file(fastaFai)
        file(cnvkit_ref)
        // file(targetBED)
    
    output:
    tuple val("cnvkit_single"), idPatient, idSample, file("*")

    when: params.cnvkit_ref && 'cnvkit_single' in tools

    script:
    
    """
    init.sh
    cnvkit.py batch ${bam} \
        --reference ${cnvkit_ref} \
        --scatter \
        --diagram
    """
}