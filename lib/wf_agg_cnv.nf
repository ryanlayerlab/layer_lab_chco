process manta_to_bed{
    tag {idPatient + "-" + idSample}
    label 'container_llab'

    //publishDir "${params.outdir}/VariantCalling/${idSample}/MantaIntermediates/", mode: params.publish_dir_mode

    input:
    tuple caller, idPatient, idSample, file(vcfgz), file(vcfgztbi)
    file(exon_file)

    output:
    tuple val('Manta'), idPatient, idSample, file("${idSample}.bed")

    when: ! ('chco_qc' in _skip_qc)  && params.bait_bed

    script:
    """
    bedtools intersect -wb -b Manta_${idSample}.tumorSV.vcf.gz -a $exon_file > temp.tsv
    manta_to_bed.py $idSample ${idSample}.bed temp.tsv
    """
}

process savvy_to_bed{
    tag {idPatient + "-" + idSample}
    label 'container_llab'

    publishDir "${params.outdir}/VariantCalling/SavvyIntermediates/", mode: params.publish_dir_mode

    input:
    file(SavvycnvResults)
    file(exon_file)

    output:
    file("*savvy.bed")

    script:
    """
    # the file is hard coded for now, change before going live to $SavvycnvResults/cnv_list.csv
    # remove chr from chromosome
    cat /scratch/Shares/CHCO/workspace/cna_positive_wes/results/savvycnv/VariantCalling/SavvycnvResults/cnv_list.csv | sed 's/chr//g' | bedtools intersect -wb -a $exon_file -b stdin > savvy-temp.tsv
    savvy_to_bed.py savvy-temp.tsv 
    ls
    """
}

process cnvkit_to_bed{
    tag {idPatient + "-" + idSample}
    label 'container_llab'

    publishDir "${params.outdir}/VariantCalling/${idSample}/ReCalToMarkedRaw/", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(cnr)
    file(exon_file)

    output:
    tuple val('cnvkit'), idPatient, idSample, file("${idSample}.bed")

    script:
    """
        # strip the first line # the next line is hard coded and needs to be fix before production as does the input parameter
        sed '1 d' /scratch/Shares/CHCO/workspace/cna_positive_wes/results/cnvkit/VariantCalling/${idSample}/CNVKit/${idSample}.cnr | bedtools intersect -wb -b stdin -a $exon_file > cnvtemp.tsv
        # sed '1 d' $cnr | bedtools intersect -wb -b stdin -a $exon_file > cnvtemp.tsv
        python cnv-kit_to_bed.py $idSample ${idSample}.bed cnvtemp.tsv
    """
}

process combine_callers{
    tag {idPatient + "-" + idSample}
    label 'container_llab'

    publishDir "${params.outdir}/VariantCalling/${idSample}/AllCNVCallers/", mode: params.publish_dir_mode

    input:
    //tuple caller, idPatient, idSample, file(manta_bed)
    tuple val(idSample)
    file(savvy_beds)
    // add input for CNVKIT eventually

    output:
    tuple idSample, file("${idSample}.multi_caller.bed")

    script:
    """
    FILE="${idSample}.coverageBinner_savvy.bed"
    if test -f "\$FILE"; then
        cat ${idSample}.coverageBinner_savvy.bed >> temp.bed
    fi 

    if test -f "temp.bed"; then
        #cnvkit_file=\$(find )

        grep -v BND temp.bed | awk '(\$2 <= \$3)' >  filtered_temp.bed
        echo "2"
        cat filtered_temp.bed | sort -k1,1V -k2,2n -k3,3n > tripple_sorted.bed
        echo "3"
        bedtools cluster -i tripple_sorted.bed > clustered_test.bed
        echo "4"
        agg_cluster.py clustered_test.bed > ${idSample}.multi_caller.bed
        echo "5"
    else
        touch ${idSample}.multi_caller.bed
    fi
    """
}

process combine_samples{
    tag {idPatient + "-" + idSample}
    label 'container_llab'

    publishDir "${params.outdir}/VariantCalling/CombineCNV/", mode: params.publish_dir_mode

    input:
    file(all_files)
    file(example_vcf)
    file(fasta)
    file(fastaFai)
    file(dict)

    output:
    file("aggregated_multi_sample_multi_caller.bed")
    file("*.vcf")
    path "cnv_all_samples.vcf", emit: cnv_all_samples_vcf

    script:
    """
    for file in ${all_files}; do
        cat \${file} >> temp.bed
    done

    #sort it and cluster it
    bedtools sort -i temp.bed > sorted_temp.bed
    cat sorted_temp.bed | sort -k1,1V -k2,2n -k3,3n > tripple_sorted.bed
    bedtools cluster -i tripple_sorted.bed > clustered_test.bed

    multi_sample_agg_cluster.py clustered_test.bed > aggregated_multi_sample_multi_caller.bed
    multi_caller_single_sample_bed_to_vcf.py --bed aggregated_multi_sample_multi_caller.bed --example_vcf $example_vcf --ref $fasta > cnv_all_samples.vcf
    """
}
