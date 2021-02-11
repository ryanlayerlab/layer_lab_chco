include {wf_gather_mapped_partial_reads} from './wf_gather_mapped_partial_reads'
include {MarkDuplicates} from './wf_mark_duplicates'
workflow wf_filter_and_gather_mapped_partial_reads{
    
    take:
       take: _partial_mapped_reads
       take: _fasta
        take: _fasta_fai
        take: _bwa_index
        take: _bam_suffix
    main:
        // bam_mapped
        // .dump('bam_mapped_input: ')
        // Run the bam QC filtering
        FilterBamRead1(_partial_mapped_reads)
        FilterBamRead2(_partial_mapped_reads)
        filtered_reads = FilterBamRead1.out.filtered_reads
                        .mix(FilterBamRead2.out.filtered_reads)
                        .groupTuple(by: [0,1,2])
                        // .dump(tag: 'filtered_reads_merge: ')
        // Merge Reads1 and Reads2
        MergeFilteredBamReads(filtered_reads)
        // Here we filter out seconday and supplemntary reads
        FilterOutSecondaryAndSupplementaryReads(
            MergeFilteredBamReads.out.filtered_bam)
        wf_gather_mapped_partial_reads(
            FilterOutSecondaryAndSupplementaryReads.out.filtered_bam_mapped,
            _fasta,
            _fasta_fai,
            _bwa_index,
            _bam_suffix) // we pass an empty string as  the bam name suffix
        
    emit:
        filtered_bams = wf_gather_mapped_partial_reads.out.bams_mapped
} // end of wf_map_reads
// STEP 1.4 FILTERING BAMS ON THE BASIS OF QUALITY SCORES AND UNIQUELY MAPPED READS. ORIGINAL COMMANDS CAME FROM REBECCA BERNARD

process FilterBamRead1 {
    label 'cpus_32'
    label 'container_llab'
    tag {idPatient + "-" + idSample + "_" + idRun}

    input:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}.bam")
        
    output:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_filtered_r1.bam"), emit: filtered_reads

    when: params.filter_bams
    script:
    if( params.remove_supplementary_reads)
        """
        init.sh
        sambamba view -t ${task.cpus} -h \
            -F "(first_of_pair and mapping_quality >=${params.bam_mapping_q} \
                and not ([XA] != null or [SA] != null)) \
                or second_of_pair" \
                "${idSample}_${idRun}.bam" \
            | samtools sort -n --threads ${task.cpus} \
            | samtools fixmate - - \
            | samtools view -h -f0x02 > "${idSample}_${idRun}_filtered_r1.bam"
        """

    else
        """
        init.sh
        sambamba view -t ${task.cpus} -h \
            -F "(first_of_pair and mapping_quality >=${params.bam_mapping_q}) \
                or second_of_pair" \
                "${idSample}_${idRun}.bam" \
            | samtools sort -n --threads ${task.cpus} \
            | samtools fixmate - - \
            | samtools view -h -f0x02 > "${idSample}_${idRun}_filtered_r1.bam"
        """
}

process FilterBamRead2 {
    label 'cpus_32'
    label 'container_llab'
    tag {idPatient + "-" + idSample + "_" + idRun}

    input:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}.bam")
        
    output:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_filtered_r2.bam"), emit: filtered_reads

    when: params.filter_bams
    script:
    if( params.remove_supplementary_reads)
        """
        init.sh
        sambamba view -t ${task.cpus} -h \
            -F "(second_of_pair and mapping_quality >=${params.bam_mapping_q} \
                and not ([XA] != null or [SA] != null)) \
                or first_of_pair" \
                "${idSample}_${idRun}.bam" \
            | samtools sort -n --threads ${task.cpus} \
            | samtools fixmate - - \
            | samtools view -h -f0x02 > "${idSample}_${idRun}_filtered_r2.bam"
        """
    else
        """
        init.sh
        sambamba view -t ${task.cpus} -h \
            -F "(second_of_pair and mapping_quality >=${params.bam_mapping_q}) \
                or first_of_pair" \
                "${idSample}_${idRun}.bam" \
            | samtools sort -n --threads ${task.cpus} \
            | samtools fixmate - - \
            | samtools view -h -f0x02 > "${idSample}_${idRun}_filtered_r2.bam"
        """
}

process MergeFilteredBamReads {
    label 'cpus_32'
    label 'container_llab'
    tag {idPatient + "-" + idSample}

    input:
        tuple idPatient, idSample, idRun, file(partial_filtered_bams)

    output:
        tuple idPatient, idSample, idRun,  file(out_bam), file("${out_bam}.bai"), emit: filtered_bam

    when: (params.filter_bams)

    script:
     out_bam = "${idSample}_${idRun}_filtered.bam"
    // bn = file(out_bam).baseName
    // out_bai = "${bn}.bai"
    """
    init.sh
    samtools merge --threads ${task.cpus} -n -c -p merged.bam ${partial_filtered_bams}
    samtools sort --threads ${task.cpus} merged.bam -o ${idSample}_${idRun}_filtered.bam
    #samtools index  --threads ${task.cpus} ${idSample}_${idRun}_filtered.bam
    samtools index   ${idSample}_${idRun}_filtered.bam
    # cleaning
    rm -f merged.bam
    """
}

process FilterOutSecondaryAndSupplementaryReads {
    label 'cpus_32'
    label 'container_llab'
    tag {idPatient + "-" + idSample}
    input:
        tuple idPatient, idSample, idRun,  file("${idSample}_${idRun}_filtered.bam"), file(bai)

    output:
        tuple idPatient, idSample, idRun, file(out_bam),  emit: filtered_bam_mapped
        // tuple idPatient, idSample, file(out_bam), file(out_bai),  emit: filtered_bam_mapped

    when: (params.filter_bams)

    script:
    out_bam = "${idSample}_${idRun}_pq${params.bam_mapping_q}.bam"
    bn = file(out_bam).baseName
    out_bai = "${bn}.bai"
    // simg = myDir = file('/scratch/Shares/layer/singularity/llab.sif')
    // simg.copyTo('/tmp/llab.sif')

    """
    #!/usr/bin/env python
    import pysam
    import shutil
    import os
    in_bam = "${idSample}_${idRun}_filtered.bam"
    # pysam.index("--threads", ${task.cpus}, in_bam)
    #pysam.index(in_bam)
    
    fd_in_bam = pysam.AlignmentFile(in_bam, "rb")
    fd_out_bam = pysam.AlignmentFile("uniq.bam", "wb", template=fd_in_bam)
    
    # Store reads to the dictionary so only unique reads stay
    reads_dict={}
    for r in fd_in_bam.fetch():
        reads_dict[(r.query_name, r.flag)] = r
    
    # iterate through the dict to write to the output bam file
    for key, read in reads_dict.items():
        fd_out_bam.write(read)
    fd_in_bam.close()
    fd_out_bam.close()
    # Sort the file by coordinates
    pysam.sort("-o", "${out_bam}", "--threads", "${task.cpus}", "uniq.bam")
    pysam.index("${out_bam}")
    shutil.move("${out_bam}.bai", "${out_bai}")
    # Cleanup
    os.unlink("uniq.bam")
    """
}