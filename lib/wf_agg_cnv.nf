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
    # the file is hard coded for now, change before going live to $SavvycnvResults/cnv_list.csv /scratch/Shares/CHCO/workspace/cna_positive_wes/results/savvycnv/VariantCalling/SavvycnvResults/cnv_list.csv
    # remove chr from chromosome
    cat $SavvycnvResults/cnv_list.csv | sed 's/chr//g' | bedtools intersect -wb -a $exon_file -b stdin > savvy-temp.tsv
    savvy_to_bed.py savvy-temp.tsv 
    ls
    """
}

process cnvkit_to_bed{
    tag {idPatient + "-" + idSample}
    label 'container_llab'

    publishDir "${params.outdir}/VariantCalling/${idSample}/ReCalToMarkedRaw/", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(cns)
    file(exon_file)

    output:
    tuple val('cnvkit'), idPatient, idSample, file("${idSample}.bed")

    script:
    """
        # strip the first line # the next line is hard coded and needs to be fix before production as does the input parameter
        sed '1 d' $cns | bedtools intersect -wb -b stdin -a $exon_file > cnvtemp.tsv
        # sed '1 d' $cns | bedtools intersect -wb -b stdin -a $exon_file > cnvtemp.tsv
        cnv-kit_to_bed.py $idSample ${idSample}.bed cnvtemp.tsv
    """
}

process combine_callers{
    tag {idPatient + "-" + idSample}
    label 'container_llab'

    publishDir "${params.outdir}/VariantCalling/${idSample}/AllCNVCallers/", mode: params.publish_dir_mode

    input:
    file(savvy_beds)
    tuple caller, idPatient, idSample, file(cnvkit_bed)
    // add input for CNVKIT eventually

    output:
    tuple idSample, file("${idSample}.multi_caller.bed"), file("${idSample}.log")

    script:
    """
    touch ${idSample}.log
    if [[ "$savvy_beds" != "not_used"  ]];
    then    
        FILE="${idSample}.coverageBinner_savvy.bed"
        if test -f "\$FILE"; then
            cat ${idSample}.coverageBinner_savvy.bed >> temp.bed
        else
            echo "# No calls from Savvy for $idSample" >> ${idSample}.log
        fi
    fi 

    if [ "$caller" != "cnvkit_not_used" ];
    then
        if [[ \$(wc -l <$cnvkit_bed) -ge 2 ]]
        then
            echo ""
        else
            echo "# No calls from CNVKit for $idSample" >> ${idSample}.log    
        fi
        echo "Adding CNVKit bed"
        cat $cnvkit_bed >> temp.bed
    else
        echo "Skipping CNVKit bed addition"
    fi

    if test -f "temp.bed"; then
        #cnvkit_file=\$(find )

        grep -v BND temp.bed | awk '(\$2 <= \$3)' >  filtered_temp.bed
        cat filtered_temp.bed | sort -k1,1V -k2,2n -k3,3n > tripple_sorted.bed
        bedtools cluster -i tripple_sorted.bed > clustered_test.bed
        agg_cluster.py clustered_test.bed > ${idSample}.multi_caller.bed
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
    file(all_logs)
    file(example_vcf)
    file(fasta)
    file(fastaFai)
    file(dict)

    output:
    file("aggregated_multi_sample_multi_caller.bed")
    file("*.vcf")
    path "cnv_all_samples.vcf", emit: cnv_all_samples_vcf
    path "cnv_all_samples.log", emit: cnv_all_samples_log

    script:
    """
    for file in ${all_files}; do
        cat \${file} >> temp.bed
    done

    for file in ${all_logs}; do
        cat \${file} >> cnv_all_samples.log
    done

    #sort it and cluster it
    bedtools sort -i temp.bed > sorted_temp.bed
    cat sorted_temp.bed | sort -k1,1V -k2,2n -k3,3n > tripple_sorted.bed
    bedtools cluster -i tripple_sorted.bed > clustered_test.bed

    multi_sample_agg_cluster.py clustered_test.bed > aggregated_multi_sample_multi_caller.bed
    multi_caller_single_sample_bed_to_vcf.py --bed aggregated_multi_sample_multi_caller.bed --example_vcf $example_vcf --ref $fasta > cnv_all_samples.vcf
    """
}
