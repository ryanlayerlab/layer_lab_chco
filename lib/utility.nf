

process IndexBamFile {
    label 'cpus_16'
    label 'container_llab'
    tag {idPatient + "-" + idSample}
    
    publishDir "${params.outdir}/Preprocessing/${idSample}/Bams/", mode: params.publish_dir_mode

    input:
        tuple idPatient, idSample, file(bam)

    output:
        tuple idPatient, idSample, file(bam), file("${bam.baseName}.bai")

    // when: !params.knownIndels

    script:
    """
    init.sh
    samtools index ${bam}
    mv ${bam}.bai ${bam.baseName}.bai
    """
}

process ConcatVCF {
    label 'cpus_8'
    label 'container_llab'

    tag {variantCaller + "-" + idSample}

    publishDir "${params.outdir}/VariantCalling/${idSample}/${output_dir}", mode: params.publish_dir_mode

    input:
        tuple variantCaller, idPatient, idSample, file(vcFiles)
        file(fastaFai)
        file(targetBED)
        val(output_file_prefix)
        val(output_file_ext)
        val(output_dir)

    output:
    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
        // tuple variantCaller, idPatient, idSample, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi"), emit: concatenated_vcf_with_index
        // tuple variantCaller, idPatient, idSample, file("*_*.vcf.gz"), emit: concatenated_vcf_without_index
        tuple variantCaller, idPatient, idSample, file("${outFile}.gz"), 
            file("${outFile}.gz.tbi"), emit: concatenated_vcf_with_index

        tuple variantCaller, idPatient, idSample, file("${outFile}.gz"), 
            emit: concatenated_vcf_without_index
            
    script:
    outFile =  "${output_file_prefix}${idSample}.${output_file_ext}"
    options = params.target_bed ? "-t ${targetBED}" : ""
    """
    init.sh
    concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outFile} ${options}
    """
}

/******************************************************************************************/
                                /* Helper functions */
/******************************************************************************************/


def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}
// Layer Lab ascii art
   def layerLabAscii() {
    // def ascii_str = 
    return
    '''
    | |                         | |         | |    
    | |     __ _ _   _  ___ _ __| |     __ _| |__ 
    | |    / _` | | | |/ _ \\ '__| |    / _` | '_ \\ 
    | |___| (_| | |_| |  __/ |  | |___| (_| | |_) |
    |______\\__,_|\\__, |\\___|_|  |______\\__,_|_.__/ 
                __/ |                            
               |___/  
    '''
    
  }

def layerLabMessage() {
  // Log colors ANSI codes
    c_reset  = params.monochrome_logs ? '' : "\033[0m";
    c_dim    = params.monochrome_logs ? '' : "\033[2m";
    c_black  = params.monochrome_logs ? '' : "\033[0;30m";
    c_red    = params.monochrome_logs ? '' : "\033[0;31m";
    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue   = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan   = params.monochrome_logs ? '' : "\033[0;36m";
    c_white  = params.monochrome_logs ? '' : "\033[0;37m";
     return """ ${c_dim}----------------------------------------------------${c_reset}
     ${c_cyan} LAYER LAB DNA Seq ANALYSIS PIPELINE ${c_reset}
     ${c_dim}----------------------------------------------------${c_reset}
     """
  
}

def helpMessage() {
  // Display help message
    log.info layerLabMessage()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --input sample.tsv -profile fiji

    Mandatory arguments:
        --input                     Path to input TSV file on mapping, recalibrate and variantcalling steps
                                    Multiple TSV files can be specified with quotes
                                    Works also with the path to a directory on mapping step with a single germline sample only
                                    Alternatively, path to VCF input file on annotate step
                                    Multiple VCF files can be specified with quotes
        -profile                    Configuration profile to use
                                    Can use multiple (comma separated)
                                    Available: conda, docker, singularity, test and more

    Options:
        --no_gvcf                    No g.vcf output from HaplotypeCaller
        --no_intervals              Disable usage of intervals
        --nucleotides_per_second      To estimate interval size
                                    Default: 1000.0
        --target_bed                Target BED file (aka primary/regions/empirical) for targeted  or whole exome sequencing
        --padded_target_bed         Target BED with extra padding that is used to calculate on_target assessment by intersecting 
                                    Recalibrated bams on these regions, and then passing through Picard CollectHsMetrics and 
                                    in conjunction to running CollectHsMetrics on Raw unmarked bams and marked recalibrated bams
        --bait_bed                  Bait BED file (aka covered/captured) for targeted or whole exome sequencing (used for GATK CollectHsMetrics)
        --step                      Specify starting step
                                    Available: Mapping, Recalibrate, VariantCalling, Annotate
                                    Default: Mapping
        --tools                     Specify tools to use for variant calling:
                                    Available: ASCAT, ControlFREEC, FreeBayes, HaplotypeCaller, DeepVariant, 
                                    Manta, mpileup, Mutect2, Mutect2_Single, gen_somatic_pon, gen_read_count_pon, 
                                    Strelka, TIDDIT
                                    and/or for annotation:
                                    snpEff, VEP, merge
                                    and for pipline validation (if you have added GIAB samples):
                                    hap_py
                                    Default: None
        --skip_qc                   Specify which QC tools to skip when running the pipeline
                                    Available: all, bamQC, BCFtools, FastQC, MultiQC, samtools, vcftools, versions
                                    Default: None
        --annotate_tools            Specify from which tools the pipeline will look for VCF files to annotate, only for step annotate
                                    Available: HaplotypeCaller, Manta, Mutect2, Strelka, TIDDIT
                                    Default: None
                                    Adds the following tools for --tools: DNAseq, DNAscope and TNscope
        --annotation_cache          Enable the use of cache for annotation, to be used with --snpEff_cache and/or --vep_cache        
        --pon_somatic               panel-of-normals VCF (bgzipped, indexed). See: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_CreateSomaticPanelOfNormals.php
        --pon_somatic_index         index of pon panel-of-normals VCF
        --pon_read_count            panel-of-normals hdf5 file. See https://gatk.broadinstitute.org/hc/en-us/articles/360040510031-CreateReadCountPanelOfNormals
        --filter_bams               Generate additional filter Bams based upon the bam_mapping_q parameter
        --bam_mapping_q             Specify the lower threshold for mapping quality (defaults to 60). All reads below this threshold will be discarded
                                    when generate the additional bams (you must specify the param --filter-bams)
        --remove_supplementary_reads When specified in addition to --filter_bams, will appy the following 
                                    sambamba expression to remove the supplementary reads 'not ([XA] != null or [SA] != null)'
    References                      If not specified in the configuration file or you wish to overwrite any of the references.
        --bwa_index                  bwa indexes
                                    If none provided, will be generated automatically from the fasta reference
        --dbsnp                     dbsnp file
        --dbsnp_index                dbsnp index
                                    If none provided, will be generated automatically if a dbsnp file is provided
        --dict                      dict from the fasta reference
                                    If none provided, will be generated automatically from the fasta reference
        --fasta                     fasta reference
        --fasta_fai                  reference index
                                    If none provided, will be generated automatically from the fasta reference
        --germline_resource          Germline Resource File
        --germline_esource_index     Germline Resource Index
                                    If none provided, will be generated automatically if a germlineResource file is provided
        --intervals                 intervals
                                    If none provided, will be generated automatically from the fasta reference
                                    Use --no_intervals to disable automatic generation
        --known_indels              knownIndels file
        --known_indels_index        knownIndels index
                                    If none provided, will be generated automatically if a knownIndels file is provided
        --snpEff_db                 snpeffDb version
        --snpEff_cache              snpeffDb cache path if you have downloaded it locally
        --species                   species for VEP
        --vep_cache                 Path to VEP cache if you have downloaded it locally
        --vep_cache_version         VEP Cache version
    Other options:
        --outdir                    The output directory where the results will be saved
        --sequencing_center         Name of sequencing center to be displayed in BAM file
        --multiqc_config            Specify a custom config file for MultiQC
        --monochrome_logs           Logs will be without colors
        --email                     Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
        --max_multiqc_email_file_size   Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
        --exome                     Specify when dealing with WES data, used when calling germline variants using Google deepvariant
        -name                       Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    AWSBatch options:
        --awsqueue                  The AWSBatch JobQueue that needs to be set when running on AWSBatch
        --awsregion                 The AWS Region for your AWS Batch job to run on
    """.stripIndent()
    // println(params.genome)
}
def printSummary(params){
    /*
    ================================================================================
                                    PRINTING SUMMARY
    ================================================================================
    */

    // Header log info
    step = params.step
    tools = params.globals.tools
    skip_qc = params.globals.skip_qc
    log.info layerLabMessage()
    def summary = [:]
    if (workflow.revision)          summary['Pipeline Release']    = workflow.revision
    summary['Run Name']          = params.custom_run_name ?: workflow.runName
    summary['Max Resources']     = "${params.max_memory} memory, ${params.max_cpus} cpus, ${params.max_time} time per job"
    if (workflow.containerEngine)   summary['Container']         = "${workflow.containerEngine} - ${workflow.container}"
    if (params.input)               summary['Input']             = params.input
    if (params.target_bed)           summary['Target BED']        = params.target_bed
    if (params.bait_bed)           summary['BAIT BED']        = params.bait_bed
    if (step)                       summary['Step']              = params.step
    if (params.tools)               summary['Tools']             = tools.join(', ')
    if (params.skip_qc)              summary['QC tools skip']     = skip_qc.join(', ')
    // if (params.filter_bams)              summary['params.filter_bams']     = params.filter_bams

    if (params.no_intervals && step != 'annotate') summary['Intervals']         = 'Do not use'
    if ('haplotypecaller' in tools)                summary['GVCF']              = params.no_gvcf ? 'No' : 'Yes'
    // if ('strelka' in tools && 'manta' in tools )   summary['Strelka BP']        = params.noStrelkaBP ? 'No' : 'Yes'
    if (params.sequencing_center)                  summary['Sequenced by']      = params.sequencing_center
    if (params.pon && 'mutect2' in tools)          summary['Panel of normals']  = params.pon
    // if (params.create_read_count_pon)          summary['Panel of normals']  = params.pon

    // summary['Save Genome Index'] = params.saveGenomeIndex ? 'Yes' : 'No'
    // summary['Nucleotides/s']     = params.nucleotidesPerSecond
    summary['Output dir']        = params.outdir
    summary['Launch dir']        = workflow.launchDir
    summary['Working dir']       = workflow.workDir
    summary['Script dir']        = workflow.projectDir
    summary['User']              = workflow.userName
    summary['genome']            = params.genome
    if (params.genomes_base)            summary['genomes base dir']   = params.genomes_base

    if (params.fasta)                 summary['fasta']                 = params.fasta
    if (params.fasta_fai)              summary['fasta_fai']              = params.fasta_fai
    if (params.fasta_gz)              summary['fasta_gz']              = params.fasta_gz
    if (params.fasta_gz_fai)          summary['fasta_gz_fai']        = params.fasta_fai
    if (params.fasta_gzi)              summary['fasta_gzi']              = params.fasta_gzi

    if (params.filter_bams)                  summary['filter_bams']                  = params.filter_bams
    if (params.bam_mapping_q)                  summary['bam_mapping_q']                  = params.bam_mapping_q
    if (params.dict)                  summary['dict']                  = params.dict
    if (params.bwa_index)              summary['bwa_index']              = params.bwa_index
    if (params.germline_resource)      summary['germline_resource']      = params.germline_resource
    if (params.germline_resource_index) summary['germline_resource_index'] = params.germline_resource_index
    if (params.intervals)             summary['intervals']             = params.intervals
    if (params.ac_loci)                summary['ac_loci']                = params.ac_loci
    if (params.ac_loci_gc)              summary['ac_loci_gc']              = params.ac_loci_gc
    if (params.chr_dir)                summary['chr_dir']                = params.chr_dir
    if (params.chr_length)             summary['chr_length']             = params.chr_length
    if (params.dbsnp)                 summary['dbsnp']                 = params.dbsnp
    if (params.dbsnp_index)            summary['dbsnp_index']            = params.dbsnp_index
    if (params.known_indels)           summary['known_indels']           = params.known_indels
    if (params.known_indels_index)      summary['known_indels_index']      = params.known_indels_index
    if (params.snpEff_db)              summary['snpEff_db']              = params.snpEff_db
    if (params.species)               summary['species']               = params.species
    if (params.vep_cache_version)       summary['vep_cache_version']       = params.vep_cache_version
    if (params.species)               summary['species']               = params.species
    if (params.snpEff_cache)          summary['snpEff_cache']          = params.snpEff_cache
    if (params.vep_cache)             summary['vep_cache']             = params.vep_cache
    if (params.cadd_InDels)            summary['cadd_InDels']          = params.cadd_InDels
    if (params.cadd_WG_SNVs)            summary['cadd_WG_SNVs']          = params.cadd_WG_SNVs
    if (params.giab_highconf_vcf)            summary['GIAB Truth Set']          = params.giab_highconf_vcf
    if (params.chco_highqual_snps)     summary['Children Colorado hig quality SNPs'] = params.chco_highqual_snps
    if (params.cadd_WG_SNVs_tbi)            summary['cadd_WG_SNVs_tbi']          = params.cadd_WG_SNVs_tbi
    

    if (workflow.profile == 'awsbatch') {
        summary['AWS Region']        = params.awsregion
        summary['AWS Queue']         = params.awsqueue
    }
    summary['Config Profile'] = workflow.profile
    if (params.config_profile_description)  summary['Config Description']  = params.config_profile_description
    if (params.config_profile_contact)      summary['Config Contact']      = params.config_profile_contact
    if (params.config_profile_url)          summary['Config URL']          = params.config_profile_url
    if (params.email) {
        summary['E-mail Address']        = params.email
        summary['MultiQC maxsize']       = params.maxMultiqcEmailFileSize
    }
    // params.properties.each { log.info "${it.key} -> ${it.value}" }
    // log.info params.dump()
    log.info summary.collect { k, v -> "${k.padRight(18)}: $v" }.join("\n")
    if (params.monochrome_logs) log.info "----------------------------------------------------"
    else log.info "\033[2m----------------------------------------------------\033[0m"
    // log.info("ByeBye")
    // println(params)
}
def sortBedIntervalsByDescendingDuration(bedIntervals){
    bedIntervals
    .map { intervalFile ->
        def duration = 0.0
        for (line in intervalFile.readLines()) {
            final fields = line.split('\t')
            if (fields.size() >= 5) duration += fields[4].toFloat()
            else {
                start = fields[1].toInteger()
                end = fields[2].toInteger()
                duration += (end - start) / params.nucleotides_per_second
            }
        }
        [duration, intervalFile]
        }.toSortedList({ a, b -> b[0] <=> a[0] })
    .flatten().collate(2)
    .map{duration, intervalFile -> intervalFile}
}

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
}

// Check parameter existence
def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    // println("passed list: $list")
    // println("real list: $realList")
    return list.every{ checkParameterExistence(it, realList) }
}

// Check if params.item exists and return params.genomes[params.genome].item otherwise
def checkParamReturnFile(item) {
    // Handle deprecation
    if (params.genomeDict && item == "dict") return file(params.genomeDict)
    if (params.genomeFile && item == "fasta") return file(params.genomeFile)
    if (params.genomeIndex && item == "fasta_fai") return file(params.genomeIndex)

    params."${item}" = params.genomes[params.genome]."${item}"
    return file(params."${item}")
}

// Define list of available tools to annotate
def defineAnnoList() {
    return [
        'HaplotypeCaller',
        'Manta',
        'Mutect2',
        'Strelka',
        'TIDDIT'
    ]
}

// Define list of skipable QC tools
def defineSkipQClist() {
    return [
        'bamqc',
        'bcftools',
        'fastqc',
        'alignment_summary',
        'insert_size_metrics',
        'markduplicates',
        'multiqc',
        'samtools',
        'hs_metrics',
        'on_target_assessment',
        'sentieon',
        'vcftools',
        'versions'
    ]
}

// Define list of available step
def defineStepList() {
    return [
        'annotate',
        'mapping',
        'markdups',
        'recalibrate',
        'variantcalling',
        'joint_genotype',
        'annotate'
    ]
}

// Define list of available tools
def defineToolList() {
    return [
        'ascat',
        'controlfreec',
        'dnascope',
        'dnaseq',
        'freebayes',
        'haplotypecaller',
        'joint_genotype',
        'deepvariant',
        // 'benchmark_dv_and_hc_against_giab',
        // 'benchmark_dv_against_hc',
        'hap_py',
        'manta',
        'merge',
        'mpileup',
        'gatk_cnv_somatic',
        'gatk_cnv_germline',
        'savvy_cnv_somatic',
        'savvy_cnv_germline',
        'cnvkit_somatic',
        'cnvkit_single',
        'mutect2',
        'mutect2_single',
        'manta_single',
        'gen_read_count_pon',
        'snpeff',
        'strelka',
        'tiddit',
        'tnscope',
        'vep'
    ]
}

// Channeling the TSV file containing BAM.
// Format is: "subject gender status sample bam bai"
// def extractBam(tsvFile) {
//     Channel.from(tsvFile)
//         .splitCsv(sep: '\t')
//         .map { row ->
//             checkNumberOfItem(row, 6)
//             def idPatient = row[0]
//             def gender    = row[1]
//             def status    = returnStatus(row[2].toInteger())
//             def idSample  = row[3]
//             def bamFile   = returnFile(row[4])
//             def baiFile   = returnFile(row[5])

//             if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
//             if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

//             return [idPatient, gender, status, idSample, bamFile, baiFile]
//         }
// }

def extractBam(tsvFile) {
    def infos = []

    def allLines = tsvFile.readLines()
    for (line in allLines){
        info("Parsing Line: ${line}")

        def trimmed = line.trim()
        def cols = trimmed.split()
        checkNumberOfItem(cols, 6)

        info("cols[0]:${cols[0]}")
        info("cols[1]:${cols[1]}")
        info("cols[2]:${cols[2]}")
        info("cols[3]:${cols[3]}")
        info("cols[4]:${cols[4]}")
        info("cols[5]:${cols[5]}")

        def idPatient  = cols[0]
        def gender     = cols[1]
        def status     = returnStatus(cols[2].toInteger())
        def idSample   = cols[3]
        def bamFile   = returnFile(cols[4])
        def baiFile   = returnFile(cols[5])
        if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
        if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

        infos.add([idPatient, gender, status, idSample, bamFile, baiFile])
    }
    return Channel.from(infos)
}

// Create a channel of germline FASTQs from a directory pattern: "my_samples/*/"
// All FASTQ files in subdirectories are collected and emitted;
// they must have _R1_ and _R2_ in their names.
def extractFastqFromDir(pattern) {
    def fastq = Channel.create()
    // a temporary channel does all the work
    Channel
        .fromPath(pattern, type: 'dir')
        .ifEmpty { error "No directories found matching pattern '${pattern}'" }
        .subscribe onNext: { sampleDir ->
            // the last name of the sampleDir is assumed to be a unique sample id
            sampleId = sampleDir.getFileName().toString()

            for (path1 in file("${sampleDir}/**_R1_*.fastq.gz")) {
                assert path1.getName().contains('_R1_')
                path2 = file(path1.toString().replace('_R1_', '_R2_'))
                if (!path2.exists()) error "Path '${path2}' not found"
                (flowcell, lane) = flowcellLaneFromFastq(path1)
                patient = sampleId
                gender = 'ZZ'  // unused
                status = 0  // normal (not tumor)
                rgId = "${flowcell}.${sampleId}.${lane}"
                result = [patient, gender, status, sampleId, rgId, path1, path2]
                fastq.bind(result)
            }
    }, onComplete: { fastq.close() }
    fastq
}

def pr_ch(ch){
    println("pr_ch() printing channel values")
    ch.view{println(it)}
}
// Extract gender and status from Channel
def extractInfos(channel) {
    def gender_map = [:]
    def status_map = [:]
    channel = channel.map{ it ->
        def idPatient = it[0]
        def gender = it[1]
        def status = it[2]
        def idSample = it[3]
        gender_map[idPatient] = gender
        status_map[idPatient, idSample] = status
        [idPatient] + it[3..-1]
    }
    [gender_map, status_map, channel]
}

def info(message){
    // log.info(message)
    if (params.debug){
        log.info(message)
    }
}

def extractFastq(tsvFile) {
    def infos = []

    def allLines = tsvFile.readLines()
    for (line in allLines){
        info("Parsing Line: ${line}")

        def trimmed = line.trim()
        def cols = trimmed.split()
        checkNumberOfItem(cols, 7)
        

        info("cols[0]:${cols[0]}")
        info("cols[1]:${cols[1]}")
        info("cols[2]:${cols[2]}")
        info("cols[3]:${cols[3]}")
        info("cols[4]:${cols[4]}")
        info("cols[5]:${cols[5]}")
        info("cols[6]:${cols[6]}")

        def idPatient  = cols[0]
        def gender     = cols[1]
        def status     = returnStatus(cols[2].toInteger())
        def idSample   = cols[3]
        def idRun      = cols[4]
        def file1      = returnFile(cols[5])
        def file2      = "null"
        file2 = returnFile(cols[6])
        if (!hasExtension(file2, "fastq.gz") && !hasExtension(file2, "fq.gz")) exit 1, "File: ${file2} has the wrong extension. See --help for more information"

        infos.add([idPatient, gender, status, idSample, idRun, file1, file2])
    }
    return Channel.from(infos)
}

def extractUnmarked(tsvFile) {
    def infos = []

    def allLines = tsvFile.readLines()
    for (line in allLines){
        info("Parsing Line: ${line}")

        def trimmed = line.trim()
        def cols = trimmed.split()
        checkNumberOfItem(cols, 6)

        info("cols[0]:${cols[0]}")
        info("cols[1]:${cols[1]}")
        info("cols[2]:${cols[2]}")
        info("cols[3]:${cols[3]}")
        info("cols[4]:${cols[4]}")
        info("cols[5]:${cols[5]}")

        def idPatient  = cols[0]
        def gender     = cols[1]
        def status     = returnStatus(cols[2].toInteger())
        def idSample   = cols[3]
        def bamFile   = returnFile(cols[4])
        if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
        def baiFile   = returnFile(cols[5])
        if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"
    
        infos.add([idPatient, gender, status, idSample, bamFile, baiFile])
    }
    return Channel.from(infos)
}


def extractDupMarked(tsvFile) {
    def infos = []

    def allLines = tsvFile.readLines()
    for (line in allLines){
        info("Parsing Line: ${line}")

        def trimmed = line.trim()
        def cols = trimmed.split()
        // checkNumberOfItem(cols, 7)
        checkNumberOfItem(cols, 6)

        info("cols[0]:${cols[0]}")
        info("cols[1]:${cols[1]}")
        info("cols[2]:${cols[2]}")
        info("cols[3]:${cols[3]}")
        info("cols[4]:${cols[4]}")
        info("cols[5]:${cols[5]}")
        // info("cols[6]:${cols[6]}")

        def idPatient  = cols[0]
        def gender     = cols[1]
        def status     = returnStatus(cols[2].toInteger())
        def idSample   = cols[3]
        def bamFile   = returnFile(cols[4])
        def baiFile   = returnFile(cols[5])
        // def recalTable = returnFile(row[6])

            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"
            // if (!hasExtension(recalTable, "recal.table")) exit 1, "File: ${recalTable} has the wrong extension. See --help for more information"            

        // infos.add([idPatient, gender, status, idSample, bamFile, baiFile, recalTable])
        infos.add([idPatient, gender, status, idSample, bamFile, baiFile])
    }
    return Channel.from(infos)
}

def extractGvcfs(tsvFile) {
    def infos = []

    def allLines = tsvFile.readLines()
    for (line in allLines){
        // info("Parsing Line: ${line}")
        def trimmed = line.trim()
        def cols = trimmed.split()
        checkNumberOfItem(cols, 6)

        def idPatient  = cols[0]
        def gender     = cols[1]
        def status     = returnStatus(cols[2].toInteger())
        def idSample   = cols[3]
        def gvcfFile   = returnFile(cols[4])
        def tbiFile   = returnFile(cols[5])
        // def recalTable = returnFile(row[6])

            if (!hasExtension(gvcfFile, "gz")) exit 1, "File: ${gvcfFile} has the wrong extension. See --help for more information"
            if (!hasExtension(tbiFile, "tbi")) exit 1, "File: ${tbiFile} has the wrong extension. See --help for more information"
            // if (!hasExtension(recalTable, "recal.table")) exit 1, "File: ${recalTable} has the wrong extension. See --help for more information"            

        // infos.add([idPatient, gender, status, idSample, bamFile, baiFile, recalTable])
        infos.add([idPatient, gender, status, idSample, gvcfFile, tbiFile])
    }
    return Channel.from(infos)
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    InputStream fileStream = new FileInputStream(path.toFile())
    InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
    BufferedReader buffered = new BufferedReader(decoder)
    def line = buffered.readLine()
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(' ')[0].split(':')
    String fcid
    int lane
    if (fields.size() == 7) {
        // CASAVA 1.8+ format
        fcid = fields[2]
        lane = fields[3].toInteger()
    } else if (fields.size() == 5) {
        fcid = fields[0]
        lane = fields[1].toInteger()
    }
    [fcid, lane]
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

// Remove .ann .gz and .vcf extension from a VCF file
def reduceVCF(file) {
    return file.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")
}

// Return status [0,1]
// 0 == Normal, 1 == Tumor
def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
}
// getvcfs to annotate
def getVCFsToAnnotate(results_dir, annotate_tools, tsv){
    vcf_to_annotate = Channel.empty()
    // vcf_no_annotate = Channel.empty()

    if (tsv_path == []) {
    // Sarek, by default, annotates all available vcfs that it can find in the VariantCalling directory
    // Excluding vcfs from FreeBayes, and g.vcf from HaplotypeCaller
    // Basically it's: results/VariantCalling/*/{HaplotypeCaller,Manta,Mutect2,SentieonDNAseq,SentieonDNAscope,SentieonTNscope,Strelka,TIDDIT}/*.vcf.gz
    // Without *SmallIndels.vcf.gz from Manta, and *.genome.vcf.gz from Strelka
    // The small snippet `vcf.minus(vcf.fileName)[-2]` catches idSample
    // This field is used to output final annotated VCFs in the correct directory
        println("annotate_tools: ${annotate_tools}")
        vcf_to_ann = 
        Channel.empty().mix(
        Channel.fromPath("${results_dir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
            .flatten().map{vcf -> ['HaplotypeCaller', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${results_dir}/VariantCalling/*/Manta/*[!candidate]SV.vcf.gz")
            .flatten().map{vcf -> ['Manta', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${results_dir}/VariantCalling/*/Mutect2/*.vcf.gz")
            .flatten().map{vcf -> ['Mutect2', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${results_dir}/VariantCalling/*/Strelka/*{somatic,variant}*.vcf.gz")
            .flatten().map{vcf -> ['Strelka', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${results_dir}/VariantCalling/*/TIDDIT/*.vcf.gz")
            .flatten().map{vcf -> ['TIDDIT', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
        )
        .filter {
            annotate_tools == [] || (annotate_tools != [] && it[0] in annotate_tools)
        }
    } else if (annotate_tools == []) {
    // Annotate user-submitted VCFs
    // If user-submitted, Sarek assume that the idSample should be assumed automatically
      vcf_to_annotate = Channel.fromPath(tsv_path)
        .map{vcf -> ['userspecified', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
    } else exit 1, "specify only tools or files to annotate, not both"

    //vcfNoAnnotate.close()
    //vcfAnnotation = vcfAnnotation.mix(vcfToAnnotate)
    return vcf_to_annotate
}
