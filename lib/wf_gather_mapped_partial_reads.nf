step = params.step
tools = params.globals.tools
status_map = params.globals.status_map
gender_map = params.globals.gender_map
// include {hasExtension} from './utility'
include {IndexBamFile} from './utility'
include {MarkDuplicates} from './wf_mark_duplicates'
workflow wf_gather_mapped_partial_reads{
    take: _partial_mapped_reads
    take: _fasta
    take: _fasta_fai
    take: _bwa_index
    take: _bam_suffix
    main:
        // _bam_suffix = ''
        // Now group the partial (mapped) bams according to the patient_id, and sample_id
        _grouped_bam_mapped = _partial_mapped_reads.groupTuple(by:[0, 1])
        _bams_for_merge = _grouped_bam_mapped
                        .map{idPatient, idSample, idRun, bams -> 
                        [ idPatient, idSample, idRun, _bam_suffix, bams ]
                        }

        // MapReads.out.bam_mapped.map{idPatient, idSample, idRun, bams -> 
        //     [ idPatient, idSample, idRun, _bam_suffix, bams ]
        // }
        // Now merge the grouped bams back to a single bam per sample
        MergeBamMapped(_bams_for_merge)
        IndexBamFile(MergeBamMapped.out)
        
        // Creating a TSV file to restart from this step
        IndexBamFile.out
        .map { idPatient, idSample, bamFile, baiFile ->
            status = status_map[idPatient, idSample]
            gender = gender_map[idPatient]
            bam = "${params.outdir}/Preprocessing/${idSample}/Bams/${idSample}${_bam_suffix}.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/Bams/${idSample}${_bam_suffix}.bai"
            bam_file = file(bam)
            bai_file = file(bai)
            "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam_file}\t${bai_file}\n"
        }.collectFile(
            name: "mapped_bam${_bam_suffix}.tsv", sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
        )
                              
    emit:
        bams_mapped = IndexBamFile.out
        // bams_mapped_qc = _out_bams_qc
} // end of wf_map_reads

/*
workflow wf_gather_mapped_reads{
    
    take:
        take: _bam_mapped
        take: _bam_suffix
    main:
        // STEP 1.5: MERGING BAM FROM MULTIPLE LANES
        (single_bams, multiple_bams) = 
        _bam_mapped.groupTuple(by:[0, 1])
        .branch{
            _: it[2].size() == 1
            __: it[2].size() > 1
        }

        // // ch_multiple_bams.subscribe{ println it}
        single_bams = single_bams.map{
            idPatient, idSample, idRun, bam ->
            [idPatient, idSample, bam]
        }
        
        // tell the mergeBamMapped to generate the .bams
        multiple_bams = 
        multiple_bams.map{idPatient, idSample, idRun, bams -> 
            [ idPatient, idSample, idRun, _bam_suffix, bams ]
        }

        MergeBamMapped(multiple_bams)
        
        _merged_bams = Channel.empty()
        _merged_bams = MergeBamMapped.out.mix(single_bams)
        IndexBamFile(_merged_bams)
        

        // Creating a TSV file to restart from this step
        _merged_bams
        .map { idPatient, idSample, bamFile ->
            status = status_map[idPatient, idSample]
            gender = gender_map[idPatient]
            bam = "${params.outdir}/Preprocessing/${idSample}/Bams/${idSample}${_bam_suffix}.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/Bams/${idSample}${_bam_suffix}.bai"
            bam_file = file(bam)
            bai_file = file(bai)
            "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam_file}\t${bai_file}\n"
        }.collectFile(
            name: "mapped_bam${_bam_suffix}.tsv", sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
        )

    emit:
        merged_bams = IndexBamFile.out
        merged_bams_qc = _merged_bams
} // end of wf_map_reads

*/
process PartitionFastQ {
    // label 'PartitionFastQ'
    label 'cpus_32'
    label 'container_llab'

    tag {idPatient + "-" + idRun}

    // publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", 
    // mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), 
        file("${idSample}_${idRun}_R2.fastq.gz")

    output:
        tuple idPatient, idSample, idRun, file("r1_split_*"), file("r2_split_*") 
        // val (tuple_list)

    when: params.split_fastq && step == 'mapping'
    
    script:
    """
    init.sh
    partition.sh \
                in=${idSample}_${idRun}_R1.fastq.gz \
                in2=${idSample}_${idRun}_R2.fastq.gz  \
                out=r1_split_%.fastq.gz \
                out2=r2_split_%.fastq.gz \
                ways=${params.split_fastq}
    """
}

// STEP 1: MAPPING READS TO REFERENCE GENOME WITH BWA MEM
process MapReads {
    label 'cpus_max'
    label 'container_llab'
    tag {idPatient + "-" + idRun}

    input:
        tuple idPatient, idSample, idRun, file(inputFile1), file(inputFile2)
        file(fasta) 
        file(fastaFai)
        file(bwaIndex) 

    output:
        // tuple idPatient, idSample, idRun, file("${idSample}_${idRun}.bam")
        // tuple idPatient, val("${idSample}_${idRun}"), file("${idSample}_${idRun}.bam")
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}.bam"), emit : bam_mapped
        // tuple idPatient, val("${idSample}_${idRun}"), file("${idSample}_${idRun}.bam"), emit : bam_mapped_BamQC
    
    // when: !(step in ['markdups','recalibrate', 'variantcalling', 'annotate'])
    script:
    // -K is an hidden option, used to fix the number of reads processed by bwa mem
    // Chunk size can affect bwa results, if not specified,
    // the number of threads can change which can give not deterministic result.
    // cf https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md
    // and https://github.com/gatk-workflows/gatk4-data-processing/blob/8ffa26ff4580df4ac3a5aa9e272a4ff6bab44ba2/processing-for-variant-discovery-gatk4.b37.wgs.inputs.json#L29
    CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
    // adjust mismatch penalty for tumor samples
    status = status_map[idPatient, idSample]
    extra = status == 1 ? "-B 3" : ""
    convertToFastq = hasExtension(inputFile1, "bam") ? "gatk --java-options -Xmx${task.memory.toGiga()}g SamToFastq --INPUT=${inputFile1} --FASTQ=/dev/stdout --INTERLEAVE=true --NON_PF=true | \\" : ""
    input = hasExtension(inputFile1, "bam") ? "-p /dev/stdin - 2> >(tee ${inputFile1}.bwa.stderr.log >&2)" : "${inputFile1} ${inputFile2}"
    """
        init.sh
        ${convertToFastq}
        bwa mem -K 100000000 -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${fasta} \
        ${input} | \
        samtools sort --threads ${task.cpus} -m 2G - > ${idSample}_${idRun}.bam
    """
}

// STEP 1.5: MERGING BAM FROM MULTIPLE LANES
process MergeBamMapped {
    label 'cpus_16'
    label 'container_llab'
    tag {idPatient + "-" + idSample}

    input:
        tuple idPatient, idSample, idRun, out_suffix, file(bams)
        // tuple idPatient, idSample, idRun, file(bam), val bam_type

    output:
        tuple idPatient, idSample,  file("${idSample}${out_suffix}.bam")

    script:
    // suffix = bams.first().minus("${idSample}_${idRun}")
    // out_file = "${idSample}${suffix}"
    """
    init.sh
    samtools merge --threads ${task.cpus} "${idSample}${out_suffix}.bam" ${bams}
    """
}
