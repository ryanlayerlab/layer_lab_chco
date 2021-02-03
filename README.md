# Layer Lab NGS Analysis Workflow (for CHCO)
 `PLEASE NOTE THAT THIS IS VERY MUCH A WORK IN PROGRESS AND THE CONFIGURATION AND THE CODE MIGHT NEED TWEAKING BEFORE RUNNING THE PIPELINE`

The pipeline is implemented in  [nextflow](https://www.nextflow.io) and implements the following workflows.
### Germline 
#### Short-Indels
1. [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)
2. [Google DeepVariant](https://github.com/google/deepvariant) 
#### Copy Number Variant Callers
1. [GATK Germline CNV Caller](https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants)

### Somatic 
#### Short-Indels
1. [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360035535892-Somatic-copy-number-variant-discovery-CNVs-)
   1. This workflow is implemnted both for tumor-normal, and tumor only (single sample mode) scenarios.

#### Copy Number Variant Callers
1. [GATK Somatic CNV Caller](https://gatk.broadinstitute.org/hc/en-us/articles/360035535892-Somatic-copy-number-variant-discovery-CNVs-)
2. [SavvyCNV](https://github.com/rdemolgen/SavvySuite) (based on off-read signals)
3. [CNVKit](https://cnvkit.readthedocs.io/en/stable/)

## Pipeline Architecture
Below are the major components of the pipeline
- Nextflow is the main glue that wraps the individual processes (individual components of analysis such as bwa-mem), and orchestrates the workflow
- Configuration files (under the `conf` directory) carries the information that nextflow uses to work out which [executor](https://www.nextflow.io/docs/latest/executor.html) (execution mechanism) to use to run the pipeline
- Individual scripts wrapped inside `nextflow processes` to carry out the actual tasks such as the *alignment* or *marking duplicate reads*, or running a *structural variant caller*
- A Singularity container is our preferred way to run the pipeline. Currently we include the singularity configuration in the subdirectory *containers*

## How to run the pipeline
In the current configuration on our Colorado University Boulder Fiji Cluster, we are using a singularity container with all the required software installed in it. So if you plan to use this, you need to have `singularity` (current version we are using is 3.1.1, but should run with others too) on your Unix *PATH* .

For a run on your machine (or any other infrastructure such as an HPC, or AWS), you need to add corresponding configuration under the `conf`. The *configuration* defines a *profile* in Nextflow lingo and needs to be passed as a parameter at the command-line when running the pipeline. At the command-prompt you will need to pass a *samples.tsv* (tab delimited) carrying at least three columns. The first column specifying the *sample name*, and the next two the paths to the *first* and the *second reads* respectively.
### Important parameters
In order to run the pipeline, the following infroamtion need to be passed at the commandline:
1. Input TSV file specifying the sample data (see an example below)
2. Reference Genome Build (GRCh37, GRCh38)
3. In case of Exom data, the target and bait bed files
4. Name(s) of the tool that needs to be run (run the pipleine with *--help* for a complete list of tools avaialbe)
5. Nextflow profile configured (see examples in the *confs* sub-directory)
   
A typical execution of the pipeline can be started as

    /scratch/Shares/layer/workspace/layer_lab_vc/main.nf \
	--name 'Testing split pipelines' \
	--input /scratch/Shares/layer/workspace/downsampled_data/cancer_center/cancer_center_two_ds.tsv \
	--target_bed /scratch/Shares/layer/ref/exome_capture_lists/S07604514_hs_hg38/S07604514_Covered.bed \
	--bait_bed /scratch/Shares/layer/ref/exome_capture_lists/S07604514_hs_hg38/S07604514_Regions.bed \
	--padded_target_bed /scratch/Shares/layer/ref/exome_capture_lists/S07604514_hs_hg38/S07604514_Padded.bed \
	--genome GRCh38 \
	-profile fiji \
	--tools mutect2_single \
	--exome \
	-resume 

#### Passing parameters to the pipeline
Nextflow (and this pipeline) allows the following three ways to pass params to a pipeline run:
1. At the *commandline*, using syntex like *--input*, or *--results-dir* (note the double dashes before param names)
2. Still at the commandline but by putting params in a *yaml* file and by passing that file as *--params my_params.yaml*. See below  an excerpt from a params.yaml file:

        input: "/Shares/layer_shared/projects/cancer_center_tiny/samples.tsv"
        runName: "cancer_center_tiny"
        sendEmail: False
        nfOptions: "-with-reports -with-trace -with-timeline -with-dag flowchart.png"
3. In various *config* files. See *conf/genomes.config* for some examples.

#### *Input.tsv*
This file contains the sample information in a seven column tsv. The columns represetns:
1. Patient ID.
2. Patient Sex (XX, XY, or ZZ for unspecified).
3. Sample type. 0 for a normal sample, and 1 for a tumor.
4. Sample Name. A single patient can have multiple samples drawn.
5. Lane Number of the flow cell carrying the smple (In future we might omit this column as the pratcie of splitting a sample in multiple lanes seems no longer relevant).
6. Read1 fastq.
7. Read2 fastq.
   
`A sample input tsv file:`

    KU1919	ZZ	0	KU1919P	1	/scratch/Shares/layer/workspace/downsampled_data/cancer_center/KU1919P.R1.fq.gz /scratch/Shares/layer/workspace/downsampled_data/cancer_center/KU1919P.R2.fq.gz
    KU1919	ZZ	1	KU1919C	1	/scratch/Shares/layer/workspace/downsampled_data/cancer_center/KU1919C.R1.fq.gz /scratch/Shares/layer/workspace/downsampled_data/cancer_center/KU1919C.R2.fq.gzL2_1.fq.gz	KU1919P/KU1919P_TD180603261_HT3LKCCXY_L2_2.fq.gz

## Porting the pipeline to your own infrastructure (HPC, AWS etc)
In order to run the pipeline on your own infrastructre, you'll need the following components:
1. All software packages used by the components that you need for your analysis. Our singularity container contains for example the GATK, BWA-MEM, Samtools, Bcftools and so on. This can be either locally installed or contained by a singularity or docker container.
2. Changes to configuration files accordingly if the software is containised (see *conf/base.config* for various directives).
3. Configuration files that provides detials about your execution envirnoment (an *executor* for nextflow, see *nextflow.config*, and  *conf/fiji.config*) for some examples.

 
## Extras
#### Building singularity container
The `containers/dna_seq` subdirectory contains the required *environment.yaml* and the singularity recipe file to generate the singularity image that we use. Our current working version of singularity is 3.1.1. Depending upon your singularity configuration you will need something like the following to generate the image:

    sudo singularity build -F layer_lab_dna_seq_gatk_4.1.4.sif layer_lab_dna_seq.def

