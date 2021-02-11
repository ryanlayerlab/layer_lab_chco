step = params.step
tools = params.globals.tools
status_map = params.globals.status_map
gender_map = params.globals.gender_map
workflow wf_mark_duplicates{
    take: _bams
    main:
        MarkDuplicates(_bams,
        'DuplicateMarked') // the output directory

         MarkDuplicates.out
        .map { idPatient, idSample, bamFile, baiFile ->
            status = status_map[idPatient, idSample]
            gender = gender_map[idPatient]
            bam = "${params.outdir}/Preprocessing/${idSample}/DuplicateMarked/${idSample}.md.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/DuplicateMarked/${idSample}.md.bai"
            bam_file = file(bam)
            bai_file = file(bai)
            "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam_file}\t${bai_file}\n"
        }.collectFile(
            name: "duplicate_marked_bams.tsv", sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
        )
    emit:
        dm_bams = MarkDuplicates.out.marked_bams
        // dm_bam_only = MarkDuplicates.out.bam_only
        // dm_bai_only = MarkDuplicates.out.bai_only
}


process MarkDuplicates {
    label 'container_llab'
    label 'cpus_max'
    tag {idPatient + "-" + idSample}

    publishDir "${params.outdir}/Preprocessing/${idSample}/${output_dir}/", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file("${idSample}.bam"), file("${idSample}.bai")
        val(output_dir)

    output:
        tuple idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bai"), emit: marked_bams
    when: step  ==  'mapping'

    script:
    """
    init.sh
    samtools sort -n --threads ${task.cpus}  -O SAM  ${idSample}.bam | \
        samblaster -M --ignoreUnmated| \
        samtools sort --threads ${task.cpus}  -O BAM > ${idSample}.md.bam

    samtools index ${idSample}.md.bam && \
        mv ${idSample}.md.bam.bai ${idSample}.md.bai
    """
}
