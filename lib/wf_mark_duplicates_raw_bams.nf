step = params.step
tools = params.globals.tools
status_map = params.globals.status_map
gender_map = params.globals.gender_map
include {MarkDuplicates} from './wf_mark_duplicates'
workflow wf_mark_duplicates_raw_bams{
    take: _bams
    main:
        MarkDuplicates(_bams,
                        'DuplicateMarkedRaw')
        MarkDuplicates.out
        .map { idPatient, idSample, bamFile, baiFile ->
            status = status_map[idPatient, idSample]
            gender = gender_map[idPatient]
            bam = "${params.outdir}/Preprocessing/${idSample}/DuplicateMarkedRaw/${idSample}.md.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/DuplicateMarkedRaw/${idSample}.md.bai"
            bam_file = file(bam)
            bai_file = file(bai)
            "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam_file}\t${bai_file}\n"
        }.collectFile(
            name: "duplicate_marked_raw_bams.tsv", sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
        )
    emit:
        dm_bams = MarkDuplicates.out.marked_bams
        // dm_bam_only = MarkDuplicates.out.bam_only
        // dm_bai_only = MarkDuplicates.out.bai_only
}