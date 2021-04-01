workflow wf_build_indexes{
    take: ch_fasta
    take: ch_dbsnp
    take: ch_germline_resource
    take: ch_known_indels
    take: ch_known_indels_index
    take: ch_somatic_pon

    main:

        BuildFastaFai(ch_fasta)
        BuildFastaGz(ch_fasta)
        BuildFastaGzFai(ch_fasta,
                        BuildFastaGz.out)
        BuildFastaGzi(BuildFastaGz.out)

        BuildBWAindexes(ch_fasta)
        BuildDict(ch_fasta)
        BuildDbsnpIndex(ch_dbsnp)
        BuildGermlineResourceIndex(ch_germline_resource)
        BuildKnownIndelsIndex(ch_known_indels)
        BuildSomaticPonIndex(ch_somatic_pon)
        fasta_fai = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : BuildFastaFai.out
        fasta_gz = params.fasta_gz ? Channel.value(file(params.fasta_gz)) : BuildFastaGz.out
        fasta_gz_fai = params.fasta_gz_fai ? Channel.value(file(params.fasta_gz_fai)) : BuildFastaGzFai.out
        fasta_gzi = params.fasta_gzi ? Channel.value(file(params.fasta_gzi)) : BuildFastaGzi.out
        bwa_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : BuildBWAindexes.out
        dict = params.dict ? Channel.value(file(params.dict)) :  BuildDict.out

        dbsnp_index = params.dbsnp ? \
                        params.dbsnp_index ? Channel.value(file(params.dbsnp_index)) \
                        : BuildDbsnpIndex.out : "null"        
        
        germline_resource_index = params.germline_resource ? \
                            params.germline_resource_index ? Channel.value(file(params.germline_resource_index)) \
                            : BuildGermlineResourceIndex.out : "null"
        
        
        known_indels_index =   params.known_indels ? \
                                params.known_indels_index ? ch_known_indels_index : \
                                BuildKnownIndelsIndex.out.collect() \
                                : "null"    
        somatic_pon_index = params.somatic_pon ? \
                                params.somatic_pon_index? ch_somatic_pon_index: \
                                BuildSomaticPonIndex.out \
                                : "null"
                                // Channel.value(file(params.somatic_pon_index)) \
                                // : BuildSomaticPonIndex.out
    emit:
        fasta_fai = fasta_fai
        fasta_gz = fasta_gz
        fasta_gz_fai = fasta_gz_fai
        fasta_gzi = fasta_gzi
        bwa_index = bwa_index
        dict = dict
        dbsnp_index = dbsnp_index
        germline_resource_index = germline_resource_index
        known_indels_index = known_indels_index
        somatic_pon_index = somatic_pon_index
        
} // end of wf_build_indices

/*
================================================================================
                                BUILDING INDEXES
================================================================================
*/

// And then initialize channels based on params or indexes that were just built
process BuildFastaGz {
    label 'container_llab'
      tag "${fasta}.gz"
    //   publishDir "$baseDir/sampleDerivatives"

      input:
      file(fasta)

      output:
      file("${fasta}.gz")
      
      when: !(params.fasta_gz)
      script:
      """
      init.sh
      bgzip -c ${fasta} > ${fasta}.gz
      """
}

process BuildFastaGzFai {
    label 'container_llab'
    tag "${fasta}.gz.fai"
    // publishDir "$baseDir/sampleDerivatives"

    input:
    file(fasta)
    file(fastagz)

    output:
    file("${fasta}.gz.fai")
    when: !(params.fasta_gz_fai)
    script:
    """
    init.sh
    samtools faidx $fastagz
    """
  }
  
process BuildFastaGzi {
    label 'container_llab'
    tag "${fasta}.gz.gzi"
    // publishDir "$baseDir/sampleDerivatives"

    input:
    file(fasta)

    output:
    file("${fasta}.gz.gzi")
    
    when: !(params.fasta_gzi)

    script:
    """
    init.sh
    bgzip -c -i ${fasta} > ${fasta}.gz
    """
  }

process BuildBWAindexes {
    label 'container_llab'
    tag {fasta}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/BWAIndex/${it}" : null }

    input:
        file(fasta)

    output:
        file("${fasta}.*")

    when: !(params.bwa_index) && params.fasta && 'mapping' in step

    script:
    """
    init.sh
    bwa index ${fasta}
    """
}



process BuildDict {
    label 'container_llab'
    tag {fasta}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/${it}" : null }

    input:
        file(fasta)

    output:
        file("${fasta.baseName}.dict")

    when: !(params.dict) && params.fasta && !('annotate' in step)

    script:
    """
    init.sh
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        CreateSequenceDictionary \
        --REFERENCE ${fasta} \
        --OUTPUT ${fasta.baseName}.dict
    """
}



process BuildFastaFai {
    label 'container_llab'
    tag {fasta}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/${it}" : null }

    input:
        file(fasta)

    output:
        file("${fasta}.fai")

    when: !(params.fasta_fai) && params.fasta && !('annotate' in step)

    script:
    """
    init.sh
    samtools faidx ${fasta}
    """
}



process BuildDbsnpIndex {
    label 'container_llab'
    tag {dbsnp}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/${it}" : null }

    input:
        file(dbsnp)

    output:
        file("${dbsnp}.tbi")

    when: !(params.dbsnp_index) && params.dbsnp && ('mapping' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools)

    script:
    """
    init.sh
    tabix -p vcf ${dbsnp}
    """
}


process BuildGermlineResourceIndex {
    label 'container_llab'
    tag {germlineResource}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/${it}" : null }

    input:
        file(germlineResource)

    output:
        file("${germlineResource}.tbi")

    when: !(params.germline_resource_index) && params.germline_resource && 'mutect2' in tools

    script:
    """
    init.sh
    tabix -p vcf ${germlineResource}
    """
}

process BuildKnownIndelsIndex {
    label 'container_llab'
    tag {knownIndels}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/${it}" : null }

    input:
        each file(knownIndels)

    output:
        file("${knownIndels}.tbi")

    when: !(params.known_indels_index) && params.known_indels && 'mapping' in step

    script:
    """
    init.sh
    tabix -p vcf ${knownIndels}
    """
}


process BuildSomaticPonIndex {
    label 'container_llab'
    tag {pon}

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_genome_index ? "reference_genome/${it}" : null }

    input:
        file(pon)

    output:
        file("${pon}.tbi")

    when: !(params.somatic_pon_index) && params.somatic_pon && ('tnscope' in tools || 'mutect2' in tools)

    script:
    """
    init.sh
    tabix -p vcf ${pon}
    """
}

