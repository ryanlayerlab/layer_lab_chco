tools = params.globals.tools

workflow wf_somalier{
    take: _dm_bam
    take: _fasta
    take: _fasta_fai
    take: _tsv_path
    take: _somalier_sites
    take: _1kg_ancestry_labels
    take: _1kg_somalier_extracted
    main:
        SomalierExtraction(_dm_bam,
                 _fasta,
                _fasta_fai,
                _somalier_sites
        )
        // Run Somalier extraction step
        _pedigree = Channel.empty()
        if ('somalier' in tools){
            if (params.pedigree){
               _pedigree = Channel.value(file(params.pedigree))
            } else{
                GenPedigreeFile(_tsv_path)
                _pedigree = GenPedigreeFile.out
            }
        }
         
        SomalierRelate(_pedigree,
                        SomalierExtraction.out.collect())
        
        SomalierAncestry(_1kg_ancestry_labels,
                        _1kg_somalier_extracted,
                        SomalierExtraction.out.collect())
    emit:
        related = SomalierRelate.out
        pedigree = _pedigree
} // end of wf_mpileup

process GenPedigreeFile {
    echo true
    label 'cpus_1'
    label 'container_llab'
    
    publishDir "${params.outdir}/Somalier/", mode: params.publish_dir_mode

    input:
        file(tsvFile)

    output:
        file("pedigree.ped")


    script:
    """
    gen_pedigree.awk ${tsvFile} > pedigree.ped
    """
}

process SomalierExtraction {
    label 'container_llab'
    label 'cpus_8'
    tag {idSample}
    
    publishDir "${params.outdir}/Somalier/extracted/", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam), file(bai)
        file(fasta)
        file(fasta_fai)
        file(somalier_sites)

    output:
        // tuple idPatient, idSample, file("${idSample}.somalier")
        file("${idSample}.somalier")
        // val("${params.outdir}/Somalier/extracted/")

    when: params.somalier_sites

    script:
    """
    init.sh
    somalier extract  --sites ${somalier_sites} -f ${fasta}  ${bam}
    """
}


process SomalierRelate {
    label 'container_llab'
    label 'cpus_16'
    // tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Somalier/relate/", mode: params.publish_dir_mode

    input:
       file(pedigree)
       file("somalier_extracted/*")
    //    file("*")

    output:
        tuple file("somalier.html"), file("somalier.pairs.tsv"), file("somalier.samples.tsv")

    when: 'somalier' in  tools

    script:
    """
    init.sh
    #mkdir somalier_extracted
    #mv *.somalier somalier_extracted/
    somalier relate --ped ${pedigree} somalier_extracted/*.somalier
    """
}

process SomalierAncestry {
    label 'container_llab'
    label 'cpus_32'
    // tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Somalier/ancestry/", mode: params.publish_dir_mode

    input:
        file(ancestry_labels_1kg_tsv)
        file(somalier_extracted_1kg)
        file("query_samples_somalier/*")

    output:
        file("somalier-ancestry.somalier-ancestry.html")
        file("somalier-ancestry.somalier-ancestry.tsv")
    
    when: 'somalier' in  tools

    script:
    """
    init.sh
    somalier ancestry --labels ${ancestry_labels_1kg_tsv} ${somalier_extracted_1kg}/*.somalier ++ query_samples_somalier/*.somalier
    """
}

