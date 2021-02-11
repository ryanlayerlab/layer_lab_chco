#!/usr/bin/env nextflow
nextflow.preview.dsl=2
// Includ helper functions from the lib/utility.nf
include {helpMessage;printSummary;
        checkNumberOfItem; checkParameterExistence;
        checkParameterList;checkParamReturnFile;
        defineAnnoList; defineSkipQClist;
        defineStepList; defineToolList;
        extractBam; extractFastqFromDir;
        extractInfos; extractFastq;
        extractUnmarked; extractDupMarked;
        flowcellLaneFromFastq; hasExtension;
        returnFile; reduceVCF;
        returnStatus; getVCFsToAnnotate;
        extractGvcfs
        } from './lib/utility' 

if (params.help) exit 0, helpMessage()
// PLATFORM = "ILLUMINA"
_THREADS = 32

// templateDir = "${workflow.projectDir}/lib/bash_templates"
// processParams()
/* Process the parameters and set the environemnt */
params.name = 'Layer Lab DNA Seq Analysis Pipeline'
params.tag = 'latest' // Default tag is latest, to be overwritten by --tag <version>
// model= params.exome ? 'WES' : 'WGS'
params.model= params.exome ? 'wes' : 'wgs'

// Check if genome exists in the config file
if (params.genomes && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the genomes.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

/* Define a global map (dictionary) to hold global values
   This map will be available to the sub-workflows and processes
   that are included from module files in the lib directory
*/
params.globals = [:]



stepList = defineStepList()
step = params.step ? params.step.toLowerCase() : ''

if (step.contains(',')) exit 1, 'You can choose only one step, see --help for more information'
if (!checkParameterExistence(step, stepList)) exit 1, "Unknown step specified: ${step} [valid steps: ${stepList}], see --help for more information"

toolList = defineToolList()
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(tools, toolList)) exit 1, 'Unknown tool(s), see --help for more information'
params.globals.tools = tools 
// Check if bm_dv_against_gatk selected, but deepvariant and haplotypecaller are not in the tools list
// if ( ('benchmark_dv_against_hc' in tools) && (!('haplotypecaller' in tools) || !('deepvariant' in tools) ) ) {
//     exit 1, "When benchmark_dv_against_hc selected, both HaplotypeCaller and DeepVariant needs to be in the tools list"
// }

skipQClist = defineSkipQClist()
skipQC = params.skip_qc ? params.skip_qc == 'all' ? skipQClist : params.skip_qc.split(',').collect{it.trim().toLowerCase()} : []
// subprocessed and workflows can access the skipQC through the params
params.globals.skip_qc = skipQC

if (!checkParameterList(skipQC, skipQClist)) exit 1, 'Unknown QC tool(s), see --help for more information'

anno_list = defineAnnoList()
// annotate_tools = params.annotate_tools ? params.annotate_tools.split(',').collect{it.trim().toLowerCase()} : []
annotate_tools = params.annotate_tools ? params.annotate_tools.split(',').collect{it.trim()} : []
params.globals.annotate_tools = annotate_tools
if (!checkParameterList(annotate_tools,anno_list)) exit 1, "Unknown tool(s) (${annotate_tools}) to annotate, see --help for more information"

// Has the run name been specified by the user?
// This has the bonus effect of catching both -name and --name
custom_run_name = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) custom_run_name = workflow.runName
params.globals.custom_run_name = custom_run_name
tsv_path = null
if (params.input && (hasExtension(params.input, "tsv") || hasExtension(params.input, "vcf") || hasExtension(params.input, "vcf.gz"))) tsv_path = params.input
if (params.input && (hasExtension(params.input, "vcf") || hasExtension(params.input, "vcf.gz"))) step = "annotate"

ch_input_sample = Channel.empty()
if (tsv_path) {
    tsvFile = file(tsv_path)
    switch (step) {
        case 'mapping': ch_input_sample = extractFastq(tsvFile); break
        case 'markdups': ch_input_sample = extractUnmarked(tsvFile); break
        case 'recalibrate': ch_input_sample = extractDupMarked(tsvFile); break
        case 'variantcalling': ch_input_sample = extractBam(tsvFile); break
        case 'joint_genotype': ch_input_sample = extractGvcfs(tsvFile); break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
}
// ch_input_sample.dump(tag: 'ch_input_sample')
// params.GLOBALS['step'] = step 

(gender_map, status_map, ch_input_sample) = extractInfos(ch_input_sample)

params.globals.gender_map= gender_map 
params.globals.status_map = status_map 


// ch_input_sample
//     .dump(tag: 'ch_input_sample after extracting info!')
/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each params in params.genomes, catch the command line first if it was defined
// params.fasta has to be the first one
params.fasta = params.genome && !('annotate' in step) ? params.genomes[params.genome].fasta ?: null : null
params.fasta_fai = params.genome && params.fasta ? params.genomes[params.genome].fasta_fai ?: null : null
params.fasta_gz = params.genome && !('annotate' in step) ? params.genomes[params.genome].fasta_gz ?: null : null
params.fasta_gz_fai = params.genome && params.fasta ? params.genomes[params.genome].fasta_gz_fai ?: null : null
params.fasta_gzi = params.genome && !('annotate' in step) ? params.genomes[params.genome].fasta_gzi ?: null : null
params.somalier_sites = params.genome ? params.genomes[params.genome].somalier_sites : null

// Annotation related params (snpEff)
params.snpEff_db = params.genome && ('annotate' in step) ? params.genomes[params.genome].snpEff_db ?: null : null
params.snpEff_cache = params.genome && ('annotate' in step) ? params.genomes[params.genome].snpEff_cache ?: null : null
params.species = params.genome && ('annotate' in step) ? params.genomes[params.genome].species ?: null : null
// VEP
params.vep_cache = params.genome && ('annotate' in step) ? params.genomes[params.genome].vep_cache ?: null : null
params.vep_cache_version = params.genome && ('annotate' in step) ? params.genomes[params.genome].vep_cache_version ?: null : null
// CADD
params.cadd_InDels = params.genome && ('annotate' in step) ? params.genomes[params.genome].cadd_InDels ?: null : null
params.cadd_InDels_tbi = params.genome && ('annotate' in step) ? params.genomes[params.genome].cadd_InDels_tbi ?: null : null
params.cadd_WG_SNVs = params.genome && ('annotate' in step) ? params.genomes[params.genome].cadd_WG_SNVs ?: null : null
params.cadd_WG_SNVs_tbi = params.genome && ('annotate' in step) ? params.genomes[params.genome].cadd_WG_SNVs_tbi ?: null : null
// CHCO related Pipeline validation using Illumina Hap.py package. GIAB truth sets
params.giab_highconf_vcf = params.genome ? params.genomes[params.genome].giab_highconf_vcf ?: null : null
params.giab_highconf_tbi = params.genome ? params.genomes[params.genome].giab_highconf_tbi ?: null : null
params.giab_highconf_regions = params.genome ? params.genomes[params.genome].giab_highconf_regions ?: null : null
params.chco_highqual_snps = params.genome ? params.genomes[params.genome].chco_highqual_snps ?: null : null


// The rest can be sorted
params.ac_loci = params.genome && 'ascat' in tools ? params.genomes[params.genome].ac_loci ?: null : null
params.ac_loci_gc = params.genome && 'ascat' in tools ? params.genomes[params.genome].ac_loci_gc ?: null : null
params.bwa_index = params.genome && params.fasta && 'mapping' in step ? params.genomes[params.genome].bwa_index ?: null : null
params.chr_dir = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chr_dir ?: null : null
params.chr_length = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chr_length ?: null : null
params.dbsnp = params.genome && \
                ('mapping' in step || 'markdups' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'joint_genotype' in tools || 'mutect2' in tools) \
                ? params.genomes[params.genome].dbsnp ?: null : null

params.dbsnp_index = params.genome && params.dbsnp ? params.genomes[params.genome].dbsnp_index ?: null : null
params.dict = params.genome && params.fasta ? params.genomes[params.genome].dict ?: null : null
params.germline_resource = params.genome ? params.genomes[params.genome].germline_resource ?: null : null
params.germline_resource_index = params.genome && params.germline_resource ? params.genomes[params.genome].germline_resource_index ?: null : null

// if user has not specified intervals on the commandline, pick them from the genomes.config file
if(! params.intervals){
    params.intervals = params.genome && !('annotate' in step) ? params.genomes[params.genome].intervals ?: null : null
}

params.known_indels = params.genome && ( 'mapping' in step || 'markdups' in step ) \
                            ? params.genomes[params.genome].known_indels ?: null : null

params.known_indels_index = params.genome && params.known_indels ? params.genomes[params.genome].known_indels_index ?: null : null


ch_somalier_sites = params.somalier_sites ? Channel.value(file(params.somalier_sites)) : "null"
ch_acLoci = params.ac_loci && 'ascat' in tools ? Channel.value(file(params.ac_loci)) : "null"
ch_acLoci_GC = params.ac_loci_GC && 'ascat' in tools ? Channel.value(file(params.ac_loci_GC)) : "null"
ch_chrDir = params.chr_dir && 'controlfreec' in tools ? Channel.value(file(params.chr_dir)) : "null"
ch_chrLength = params.chr_length && 'controlfreec' in tools ? Channel.value(file(params.chr_length)) : "null"
dbsnp_check = ('mapping' in step || 'joint_genotype' in step || 'controlfreec' in tools || 'haplotypecaller' in tools|| 'joint_genotype' in tools || 'mutect2' in tools)
ch_dbsnp = params.dbsnp &&  dbsnp_check ? Channel.value(file(params.dbsnp)) : "null"
// ch_dbsnp_index = params.dbsnp_index && ('mapping' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools) ? Channel.value(file(params.dbsnp_index)) : "null"
ch_fasta = params.fasta && !('annotate' in step) ? Channel.value(file(params.fasta)) : "null"
// ch_dict = params.dict ? Channel.value(file(params.dict)) : "null"
// ch_fasta_fai = params.fasta_fai && !('annotate' in step) ? Channel.value(file(params.fasta_fai)) : "null"
// ch_germline_resource = params.germline_resource && 'mutect2' in tools ? Channel.value(file(params.germline_resource)) : "null"
// ch_intervals = params.intervals && !params.no_intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : "null"
ch_germline_resource = params.germline_resource &&
                      ('mutect2' in tools || 'mutect2_single' in tools || 'gen_somatic_pon' in tools ) 
                      ? Channel.value(file(params.germline_resource)) : "null"
// if (ch_germline_resource == 'null') error "ch_germline_resource is empty!${params.germline_resource}"
ch_intervals = params.intervals && !params.no_intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : "null"

ch_read_count_pon = params.read_count_pon ? Channel.value(file(params.read_count_pon)) : "null"
ch_somatic_pon = params.somatic_pon ? Channel.value(file(params.somatic_pon)) : "null"
ch_somatic_pon_index = params.somatic_pon_index ? Channel.value(file(params.somatic_pon_index)) : "null"
ch_target_bed = params.target_bed ? Channel.value(file(params.target_bed)) : "null"
// padded target, if specified also generate CollectHsMetrics for the recal bams intersected with the padded bed
ch_padded_target_bed = params.padded_target_bed ? Channel.value(file(params.padded_target_bed)) : "null"
ch_bait_bed = params.bait_bed ? Channel.value(file(params.bait_bed)) : "null"
// knownIndels is currently a list of file for smallGRCh37, so transform it in a channel
li_known_indels = []
if (params.known_indels && ('mapping' in step || 'markdups' in step)) params.known_indels.each { li_known_indels.add(file(it)) }

li_known_indels_index = []
if (params.known_indels_index && ('mapping' in step || 'markdups' in step)) params.known_indels_index.each { li_known_indels_index.add(file(it)) }
// ch_known_indels = Channel.empty()
// ch_known_indels_index = Channel.empty()
ch_known_indels = params.known_indels && params.genome == 'smallGRCh37' ? Channel.value(li_known_indels.collect()) 
    : params.known_indels ? Channel.value(file(params.known_indels)) : "null"

ch_known_indels_index = params.known_indels_index && params.genome == 'smallGRCh37' \
    ? Channel.value(li_known_indels_index.collect()) \
    : \
    params.known_indels_index ? Channel.value(file(params.known_indels_index)) : "null"

ch_snpEff_cache = params.snpEff_cache ? Channel.value(file(params.snpEff_cache)) : "null"
// ch_snpeff_cache = params.snpeff_cache ? Channel.fromPath(params.snpeff_cache) : "null"
ch_snpEff_db = params.snpEff_db ? Channel.value(params.snpEff_db) : "null"


ch_vep_cache_version = params.vep_cache_version ? Channel.value(params.vep_cache_version) : "null"
ch_vep_cache = params.vep_cache ? Channel.value(file(params.vep_cache)) : "null"
// Optional files, not defined within the params.genomes[params.genome] scope
ch_cadd_InDels = params.cadd_InDels ? Channel.value(file(params.cadd_InDels)) : "null"
ch_cadd_InDels_tbi = params.cadd_InDels_tbi ? Channel.value(file(params.cadd_InDels_tbi)) : "null"
ch_cadd_WG_SNVs = params.cadd_WG_SNVs ? Channel.value(file(params.cadd_WG_SNVs)) : "null"
ch_cadd_WG_SNVs_tbi = params.cadd_WG_SNVs_tbi ? Channel.value(file(params.cadd_WG_SNVs_tbi)) : "null"
ch_cnvkit_ref = params.cnvkit_ref ? Channel.value(file(params.cnvkit_ref)) : "null"

// Optional CHCO files for calculating TP, FP, TN etc against the GIAB
ch_giab_highconf_vcf = params.giab_highconf_vcf ? Channel.value(file(params.giab_highconf_vcf)) : "null"
ch_giab_highconf_tbi = params.giab_highconf_tbi ? Channel.value(file(params.giab_highconf_tbi)) : "null"
ch_giab_highconf_regions = params.giab_highconf_regions ? Channel.value(file(params.giab_highconf_regions)) : "null"
ch_chco_highqual_snps = params.chco_highqual_snps ? Channel.value(file(params.chco_highqual_snps)) : "null"

// Parameters needed for QC
qc_finger_print_sites = params.finger_print_sites ? Channel.value(file(params.finger_print_sites)) : "null"
qc_extra_finger_print_sites = params.extra_finger_print_sites ? Channel.value(file(params.extra_finger_print_sites)) : "null"

/* Create channels for various indices. These channels are either filled by the user parameters or 
form inside the build_indices workflow */
ch_fasta_fai = ch_fasta_gz = ch_fasta_gzi = ch_fasta_gz_fai \
= ch_bwa_index = ch_dict = ch_dbsnp_index = ch_germline_resource_index \
=  ch_somatic_pon_index =  Channel.empty()
printSummary(params)


include {wf_get_software_versions} from './lib/wf_get_software_versions' 
include {wf_build_indexes} from './lib/wf_build_indexes' 
include {wf_build_intervals} from './lib/wf_build_intervals' 
include {wf_fastqc_fq} from './lib/wf_fastqc_fq' 
include {wf_map_partial_reads} from './lib/wf_map_partial_reads' 
include {wf_gather_mapped_partial_reads} from './lib/wf_gather_mapped_partial_reads' 
include {wf_filter_and_gather_mapped_partial_reads} from './lib/wf_filter_and_gather_mapped_partial_reads' 
include {wf_qc_bam_mapped} from './lib/wf_qc_bam_mapped' 
include {wf_mark_duplicates} from './lib/wf_mark_duplicates' 
include {wf_mark_duplicates_raw_bams} from './lib/wf_mark_duplicates_raw_bams.nf' 
include {wf_somalier_extraction} from './lib/wf_somalier_extraction' 
include {wf_recal_bam} from './lib/wf_recal_bam' 
include {wf_qc_bam_recal} from './lib/wf_qc_bam_recal' 
include {wf_deepvariant} from './lib/wf_deepvariant' 
include {wf_mpileup} from './lib/wf_mpileup' 
include {wf_mutect2_single} from './lib/wf_mutect2_single' 
include {wf_mutect2_TN} from './lib/wf_mutect2_TN' 
include {wf_haplotypecaller} from './lib/wf_haplotypecaller' 
include {wf_individually_genotype_gvcf} from './lib/wf_individually_genotype_gvcf' 
include {wf_jointly_genotype_gvcf} from './lib/wf_jointly_genotype_gvcf' 
include {wf_gatk_cnv_somatic} from './lib/wf_gatk_cnv_somatic' 
include {wf_savvy_cnv_somatic} from './lib/wf_savvy_cnv_somatic' 
include {wf_cnvkit_somatic} from './lib/wf_cnvkit_somatic' 
include {wf_manta_single} from './lib/wf_manta_single' 
include {wf_vcf_stats} from './lib/wf_vcf_stats' 
include {wf_multiqc} from './lib/wf_multiqc' 
include {ConcatVCF} from './lib/wf_haplotypecaller'
include {wf_alamut} from './lib/wf_alamut'
include {exonCoverage; onTarget; wf_raw_bam_exonCoverage; insertSize; dnaFingerprint; collectQC; wf_qc_fingerprinting_sites} from './lib/wf_quality_control'

workflow{

    wf_get_software_versions()
    wf_build_indexes(ch_fasta,
                     ch_dbsnp,
                    ch_germline_resource,
                    ch_known_indels,
                    ch_known_indels_index,
                    ch_somatic_pon)

    ch_fasta_fai =  wf_build_indexes.out.fasta_fai
    // The following three mainly used by the DeepVariant
    ch_fasta_gz =  wf_build_indexes.out.fasta_gz
    ch_fasta_gz_fai =  wf_build_indexes.out.fasta_gz_fai
    ch_fasta_gzi =  wf_build_indexes.out.fasta_gzi
    ch_bwa_index =  wf_build_indexes.out.bwa_index
    ch_dict =  wf_build_indexes.out.dict
    ch_dbsnp_index = wf_build_indexes.out.dbsnp_index
    ch_germline_resource_index = wf_build_indexes.out.germline_resource_index
    ch_known_indels_index = wf_build_indexes.out.known_indels_index.collect()
    ch_somatic_pon_index = wf_build_indexes.out.somatic_pon_index
    
    wf_build_intervals(ch_fasta_fai)
    ch_bed_intervals = wf_build_intervals.out.bed_intervals
    if (params.no_intervals && step != 'annotate') 
        ch_bed_intervals = Channel.from(file("no_intervals.bed"))
    
    // FastQCFQ(ch_input_sample)     
    ch_bam_mapped = Channel.empty()
    ch_bam_mapped_raw = Channel.empty()
    ch_fastqc_report = Channel.empty()
    
    if (step == 'mapping'){
        wf_fastqc_fq(ch_input_sample)
        // split the fastqs and map the partial fastqs
        wf_map_partial_reads(ch_input_sample,
                    ch_fasta,
                    ch_fasta_fai,
                    ch_bwa_index
        )
        ch_partial_mapped_reads = wf_map_partial_reads.out
        // gather the partial mapped reads to merge them per sample
        wf_gather_mapped_partial_reads(ch_partial_mapped_reads,
                                ch_fasta,
                                ch_fasta_fai,
                                ch_bwa_index,
                                '')
        ch_bam_mapped = wf_gather_mapped_partial_reads.out.bams_mapped
        // Now we check if we need to filter the bams
        if (params.filter_bams){
            // In this case we'll use the filtered bams for most downstream analyses
            // We still keep the raw (unfiltered bams) around for copy number callers
            // ch_partial_mapped_reads.dump(tag:'target for filtering:')
            wf_filter_and_gather_mapped_partial_reads(ch_partial_mapped_reads,
                                ch_fasta,
                                ch_fasta_fai,
                                ch_bwa_index,
                                "_pq${params.bam_mapping_q}")
            // Since ch_bam_mapped (the raw bams) wont be used often for downstream analysis
            // we keep them aside but only after duplcate marking
            wf_mark_duplicates_raw_bams(ch_bam_mapped)
            
            ch_bam_mapped = wf_filter_and_gather_mapped_partial_reads.out
        }
            
        ch_fastqc_report = wf_fastqc_fq.out.fastqc_reports.collect()
    } 

    // If the pipeline is being started from the step 'markdups'
    // use bams from the provided tsv
    if (step == 'markdups'){        

       ch_bam_mapped = ch_input_sample
    }

    //QC raw bams
    wf_qc_bam_mapped(ch_bam_mapped,
                        ch_target_bed)

    ch_bam_marked = Channel.empty()
    
    ch_bam_marked = Channel.empty()
    if (!(step in ['recalibrate', 'variantcalling', 'joint_genotype', 'annotate'])){
            wf_mark_duplicates(ch_bam_mapped)
            ch_bam_marked = wf_mark_duplicates.out.dm_bams
    }
    // Run Somalier extraction step
    wf_somalier_extraction(ch_bam_marked,
                            ch_fasta,
                            ch_fasta_fai,
                            ch_somalier_sites)
    // GATK Base Quality Score Recalibration
    wf_recal_bam(
            ch_bam_marked, // recalibrated bams
            ch_bed_intervals,
            ch_fasta,
            ch_fasta_fai,         
            ch_dict,
            ch_dbsnp,
            ch_dbsnp_index,
            ch_known_indels,
            ch_known_indels_index
        )
    ch_bam_recal = wf_recal_bam.out.bams_recal
    // Run QC metrics generation processes on the recalibrated bams

    wf_qc_bam_recal(
            ch_bam_recal,
            // ch_bam_recal_on_target,
            ch_target_bed,
            ch_bait_bed,
            ch_fasta,
            ch_fasta_fai,         
            ch_dict
            )

    // QC stats that go into final QC excel report
    exonCoverage(ch_bam_recal,ch_fasta,ch_fasta_fai,ch_dict,ch_target_bed,ch_bait_bed,"recal")
    onTarget(ch_bam_recal,ch_fasta,ch_fasta_fai,ch_dict,ch_target_bed,ch_padded_target_bed)
    wf_raw_bam_exonCoverage(wf_mark_duplicates_raw_bams.out.dm_bams,ch_fasta,ch_fasta_fai,ch_dict,ch_target_bed,ch_bait_bed)
    insertSize(ch_bam_recal)
    wf_qc_fingerprinting_sites(ch_bam_recal,qc_extra_finger_print_sites)
    dnaFingerprint(ch_bam_recal,qc_finger_print_sites,"Normal")

/* At this point we have the following bams:
a) raw unmarked bams
b) marked bams)
c) recalibrated bams
*/
    // Use Duplicated Marked Bams for DeepVariant
    wf_deepvariant(
                ch_bam_marked,
                ch_target_bed,
                ch_fasta,
                ch_fasta_fai,
                ch_fasta_gz,
                ch_fasta_gz_fai,
                ch_fasta_gzi
                )
    // Create a cartesian product of the bams and the interval files
    // for parallelization
    ch_int_bam_marked = ch_bam_marked.combine(ch_bed_intervals)
    
    wf_mpileup(
        ch_int_bam_marked,
        ch_fasta,
        ch_fasta_fai
    )

    ch_bam_for_vc = Channel.empty()
    if ((step == 'variantcalling') || (step ==  'joint_genotype') ){
       ch_bam_for_vc = ch_input_sample 
    } else{
        ch_bam_for_vc = ch_bam_recal
    }

    // Handle the Mutect2 related workflows

    
    
    // We need to prepare samples going to the mutect2 single sample mode.
    // This depends upon what is included in the params.tools. 
    // Is it just the mutect2_single
    // or gen_somatic_pon or both
    // For generating a somatic pon, we only need normal bams.
    // When running mutect2_single without the intent of generating a somatic_pon,
    // run it on all bams (tumors plus normals)
    
    // separate BAM by status, i.e normal  or tumor
    // bamNormal = Channel.empty()
    // bamTumor = Channel.empty ()
  

    (ch_bam_normal, ch_bam_tumor) = 
        ch_bam_for_vc.branch{
            _:  status_map[it[0], it[1]] == 0
            __: status_map[it[0], it[1]] == 1
        }

    // A channel for holding the cartesian product of interval files and //
    // bams for mutect2 in single sample mode
    ch_int_bam_mutec2_single = Channel.empty()
    if ('mutect2_single' in tools ) // use all bams (normal plus tumor)
        ch_int_bam_mutec2_single = ch_bam_for_vc.combine(ch_bed_intervals)
    // else if ('gen_somatic_pon' in tools && !('mutect2_single' in tools)){
    //     // use normal only
    //     ch_int_mutec2_single_bams = bam_normal.combine(ch_bed_intervals)
    // }
    
    //  ch_int_mutec2_single.dump(tag: 'ch_int_mutec2_single_bams: ')
    
    wf_mutect2_single(
        ch_int_bam_mutec2_single,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_index,
        ch_target_bed
    )

    
    // bam_normal.dump(tag: 'bam_normal: ')
    // bam_tumor.dump(tag: 'bam_tumor: ')

    // Crossing Normal and Tumor to get a T/N pair for Somatic Variant Calling
    // Remapping channel to remove common key idPatient
    // normal[0], and tumor[0] both are idPatient, so we keep just one at the start of the tuple
    ch_bam_pair = ch_bam_normal.cross(ch_bam_tumor).map {
        normal, tumor ->
        [normal[0], normal[1], normal[2], normal[3], tumor[1], tumor[2], tumor[3]]
    }
    // pair_bam.dump(tag:'BAM Somatic Pair')
    // take a cross product of the pair_bams with intervals for parrelization
    ch_int_bam_pair = ch_bam_pair.combine(ch_bed_intervals)
    // ch_int_pair_bam.dump(tag:'ch_int_pair_bam')
   // Mutect2 Tumor Normal workflow
    wf_mutect2_TN(
        ch_int_bam_pair,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_index,
        ch_somatic_pon,
        ch_somatic_pon_index,
        ch_target_bed
    )

    ch_int_bam_recal = ch_bam_for_vc
                    .combine(ch_bed_intervals)
  
    wf_haplotypecaller(
        ch_int_bam_recal,
        ch_fasta,
        ch_fasta_fai,         
        ch_dict,
        ch_target_bed,
        ch_dbsnp,
        ch_dbsnp_index
    )
   

    // individually genotype gvcfs
    wf_individually_genotype_gvcf(
        wf_haplotypecaller.out.gvcf_GenotypeGVCFs,
        ch_fasta,
        ch_fasta_fai,
        ch_dict,
        ch_dbsnp,
        ch_dbsnp_index,
        ch_target_bed
    )
// Check per sample gvcfs and prepare an input for the joint genotyping step
// If 'haplotypecaller' was specified as one of the tools, get its output and use that
// as an input; otherwise extract gvcfs from the input samples.tsv
    ch_gvcfs_for_joint_genotyping = Channel.empty()
    if ('haplotypecaller' in tools && 'joint_genotype' in tools){
        ch_gvcfs_for_joint_genotyping = wf_haplotypecaller.out.gvcf_per_sample
    }else if (step == 'joint_genotype' && 'joint_genotype' in tools) {
        ch_gvcfs_for_joint_genotyping = ch_input_sample
                                        .map{idP, idS, gvcf, tbi ->
                                            ['HaplotypeCallerGVCF', idP,idS,gvcf,tbi]
                                            }
    }

    // ch_gvcfs_for_joint_genotyping.dump(tag: 'ch_gvcfs_for_joint_genotyping')
    // Now gather all gvcfs for input to jointly_genotype_gvcf
    ch_sample_ids = ch_gvcfs_for_joint_genotyping
                    .map{vc, idP, idS, gvcf,tbi ->
                        [idS]
                        }
                        
    ch_gvcfs = ch_gvcfs_for_joint_genotyping.map{vc, idP, idS, gvcf,tbi ->
                        gvcf
                        }
                        .collect()
    
    ch_tbis = ch_gvcfs_for_joint_genotyping.map{vc, idP, idS, gvcf,tbi ->
                        tbi
                        }
                        .collect()
    
    wf_jointly_genotype_gvcf(
        ch_sample_ids,
        ch_gvcfs,
        ch_tbis,       
        ch_fasta,
        ch_fasta_fai,         
        ch_bed_intervals,
        ch_target_bed,
        ch_dict,
        ch_dbsnp,
        ch_dbsnp_index
    )

    wf_gatk_cnv_somatic (
        // wf_mark_duplicates.out.dm_bams,
        ch_bam_for_vc,
        ch_target_bed,
        ch_fasta,
        ch_fasta_fai,         
        ch_dict,
        ch_read_count_pon
    )
    
    wf_savvy_cnv_somatic(ch_bam_for_vc)
    
    (ch_normal_md_bam, ch_tumor_md_bam) = 
        ch_bam_marked.branch{
            _:  status_map[it[0], it[1]] == 0
            __: status_map[it[0], it[1]] == 1
        }
    ch_normal_md_bam_only = ch_normal_md_bam.map{
                            idP,idS,bam,bai -> 
                            bam
                            }.collect()
    ch_normal_md_bai_only = ch_normal_md_bam.map{
                            idP,idS,bam,bai -> 
                            bai
                            }.collect()
    ch_tumor_md_bam_only = ch_tumor_md_bam.map{
                            idP,idS,bam,bai -> 
                            bam
                            }.collect()
    ch_tumor_md_bai_only = ch_tumor_md_bam.map{
                            idP,idS,bam,bai -> 
                            bai
                            }.collect()
    
    wf_cnvkit_somatic( ch_normal_md_bam_only,
                        ch_normal_md_bai_only,
                        ch_tumor_md_bam_only,
                        ch_tumor_md_bai_only,
                        ch_fasta,
                        ch_fasta_fai,
                        ch_target_bed
                        )
    wf_manta_single(ch_bam_for_vc,
                    ch_target_bed,
                    ch_fasta,
                    ch_fasta_fai
                    )
    // wf_cnvkit_cnv(
    //     ch_bams.collect(),
    //     ch_fasta,
    //     ch_target_bed
    //     )
    
    // wf_tiddit_gm_sv(
    //     ch_bams,
    //     ch_fasta,
    //     ch_fasta_fai
    //     )

    //wf_vcf_stats(wf_deepvariant.out.vcf,
    //    wf_jointly_genotype_gvcf.out.vcfs_with_indexes
    //    )
    wf_vcf_stats(wf_deepvariant.out.vcf,
        wf_haplotypecaller.out.gvcf_per_sample
       )

    wf_multiqc(
        wf_get_software_versions.out,
        ch_fastqc_report,
        wf_qc_bam_mapped.out.bam_qc,
        wf_qc_bam_recal.out.samtools_stats,
        wf_qc_bam_recal.out.alignment_summary_metrics,
        wf_qc_bam_recal.out.insert_size_metrics, 
        wf_qc_bam_recal.out.hs_metrics.ifEmpty([]),
        wf_vcf_stats.out.bcfootls_stats,
        wf_vcf_stats.out.vcfootls_stats
    )
    wf_alamut(wf_haplotypecaller.out.gvcf_per_sample)
    //wf_alamut(wf_jointly_genotype_gvcf.out.vcfs_with_indexes)
    // wf_alamut(wf_jointly_genotype_gvcf.out.vcfs_with_indexes)
    collectQC(file(tsv_path), params.outdir,exonCoverage.out.files.collect(),wf_raw_bam_exonCoverage.out.raw_onTarget.collect(),insertSize.out,dnaFingerprint.out.collect(),wf_vcf_stats.out.bcfootls_stats.collect(),wf_alamut.out)
} // end of workflow



workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
