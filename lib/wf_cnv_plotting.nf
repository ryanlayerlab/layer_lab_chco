workflow wf_cnv_data_prepossessing{
    take: _bams_unfiltered
    take: _probes
    take: _ref
    take: _ref_fai
    take: _ref_dict
    take: _genes_file

    main:
    mosdepth(_bams_unfiltered)
    gzip_probes(_probes)
    count_reads(_bams_unfiltered)
    get_total_reads(count_reads.out.map{idP, idS, counts -> [counts]}.collect())
    get_probe_reads_per_million(mosdepth.out,get_total_reads.out,gzip_probes.out)
    exon_coverage_rates(mosdepth.out, gzip_probes.out)
    get_regions_zscores(exon_coverage_rates.out.map{idP, idS, bed_gz, bed_gz_tbi -> [bed_gz, bed_gz_tbi]}.collect())
    get_adj_zscore(exon_coverage_rates.out, get_regions_zscores.out)
    get_probe_cover_mean_std(exon_coverage_rates.out.map{idP, idS, bed_gz, bed_gz_tbi -> [bed_gz, bed_gz_tbi]}.collect())
    merge_adj_scores(get_adj_zscore.out.map{idP, idS, bed_gz, bed_gz_tbi -> [bed_gz, bed_gz_tbi]}.collect())

    println(get_adj_zscore.out.map{idP, idS, bed_gz, bed_gz_tbi -> [bed_gz, bed_gz_tbi]}.collect().view())

    // allele balance portion
    collect_allele_counts(_bams_unfiltered,_probes,_ref,_ref_fai,_ref_dict)
    agg_allele_counts(collect_allele_counts.out, _probes)
    merge_all_allele_counts(agg_allele_counts.out.map{ idP, idS, counts -> [counts]}.collect(), _probes)
    label_exons(_genes_file)

    emit:
    adj_probe_scores = merge_adj_scores.out
    allele_balance = merge_all_allele_counts.out
    labeled_exons = label_exons.out
    probe_cover_mean_std = get_probe_cover_mean_std.out
    
} //  end of wf_cnv_coverage_depth

process mosdepth {
    // mosdepth needs to be added to the llab container
    //label 'container_llab'
    label 'cpus_1'

    // tag {idSample}
    publishDir "${params.outdir}/CNV_Plotting/${idSample}/Mosdepth", mode: params.publish_dir_mode
    input:
        tuple idPatient, idSample, file(bam), file(bai)

    output:
        tuple idPatient, idSample, file("${idSample}.per-base.bed.gz"), file("${idSample}.per-base.bed.gz.tbi")


    // when: 'haplotypecaller' in tools

    script:
    """
    mosdepth $idSample $bam
    tabix -p bed ${idSample}.per-base.bed.gz
    """
}

// Get number of reads per bam
process count_reads {
    label 'container_llab'
    label 'cpus_1'
    publishDir "${params.outdir}/CNV_Plotting/${idSample}/ReadCounts", mode: params.publish_dir_mode 
    
    input:
        tuple idPatient, idSample, file(bam), file(bai)
 
    output:
        tuple idPatient, idSample, file("${idSample}.num_reads.txt")
    
    script:
    """
    samtools view -c -F 260 $bam > ${idSample}.num_reads.txt
    """

}


process get_probe_reads_per_million{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/${idSample}/RPM", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(per_base_bed_gz), file(per_base_bed_gz_tbi)
    file(total_reads) // out put from get_total_reads
    tuple file(probes_gz), file(probes_gz_tbi)

    output:
    tuple idPatient, idSample, file("${idSample}.probe.rpm_rate.bed.gz")  // TBD each output needs to be named for each sample

    script:
    """
    get_rpm_rates.py \
    -m $per_base_bed_gz \
    --regions_file $probes_gz \
    --num_reads $total_reads \
    | bgzip -c > ${idSample}.probe.rpm_rate.bed.gz
    """
}

process get_total_reads{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/TotalReads", mode: params.publish_dir_mode

    input:
    file(num_read_files) // all the out puts from get_probe_reads_per_million

    output:
        file("total_read.txt")

    script:
    """
    cat *.num_reads.txt > total_read.txt
    """
}

process gzip_probes{
    label 'container_llab'
    label 'cpus_1'

    input:
    file(probes)
    
    output:
    tuple file("${probes}.sorted.gz"), file("${probes}.sorted.gz.tbi")
    
    """
    bedtools sort -i $probes > ${probes}.sorted
    bgzip ${probes}.sorted
    tabix ${probes}.sorted.gz -p bed
    """
}

process exon_coverage_rates {
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/${idSample}/ExonCoverageRates", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(per_base_bed_gz), file(per_base_bed_gz_tbi) // mosdepth output
    tuple file(probes_gz), file(probes_gz_tbi) // gzip_probes output 

    output:
    tuple idPatient, idSample, file("${idSample}.probe.coverage_rate.bed.gz"), file("${idSample}.probe.coverage_rate.bed.gz.tbi")
        
    script:
    """
    get_coverage_rates.py \
    -m $per_base_bed_gz \
    --exons_file $probes_gz \
    | bgzip -c > ${idSample}.probe.coverage_rate.bed.gz
    tabix -p bed ${idSample}.probe.coverage_rate.bed.gz
    """
}

process get_regions_zscores {

    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/RegionZscores", mode: params.publish_dir_mode

    input:
    file(probe_coverage_rate_bed_gz) // all the outputs from exon_coverage_rates

    output:
    tuple file("probe.cover.mean.stdev.bed.gz"), file("probe.cover.mean.stdev.bed.gz.tbi"), file("probe.cover.mean.stdev.bed")

    script:
    """
    get_regions_zscores.py -r ./ | bgzip -c > probe.cover.mean.stdev.bed.gz
    cp probe.cover.mean.stdev.bed.gz copy_probe.cover.mean.stdev.bed.gz
    gunzip copy_probe.cover.mean.stdev.bed.gz 
    cp copy_probe.cover.mean.stdev.bed probe.cover.mean.stdev.bed
    tabix -p bed probe.cover.mean.stdev.bed.gz
    """   

}

process get_probe_cover_mean_std {

    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/ProbeCoverMeanSTD", mode: params.publish_dir_mode

    input:
    tuple file(probe_coverage_rate_bed_gz), file(probe_coverage_rate_bed_gz_tbi)// all the outputs from exon_coverage_rates

    output:
    tuple file("probe.cover.mean.stdev.bed.gz"), file("probe.cover.mean.stdev.bed.gz.tbi"), file("probe.cover.mean.stdev.bed")

    script:
    """
    get_regions_zscores.py -r ./ | bgzip -c > probe.cover.mean.stdev.bed.gz
    cp probe.cover.mean.stdev.bed.gz copy_probe.cover.mean.stdev.bed.gz
    gunzip copy_probe.cover.mean.stdev.bed.gz 
    cp copy_probe.cover.mean.stdev.bed probe.cover.mean.stdev.bed
    tabix -p bed probe.cover.mean.stdev.bed.gz
    """

}

process get_adj_zscore{

    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/${idSample}/AdjZscore", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(coverage_rate_bed_gz), file(coverage_rate_bed_gz_tbi) //out out from exon_coverage_rates
    tuple file(exon_coverage_bed_gz), file(exon_coverage_bed_gz_tbi) //out put from get_regions_zscores

    output:
    tuple idPatient, idSample, file("${idSample}.adj_z.bed.gz"), file("${idSample}.adj_z.bed.gz.tbi")

    script:
    """
    get_coverage_zscores.py \
        -r $coverage_rate_bed_gz \
        -s $exon_coverage_bed_gz \
        | bgzip -c > ${idSample}.adj_z.bed.gz
    tabix -p bed ${idSample}.adj_z.bed.gz
    """
}

process merge_adj_scores{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/MergeAdjZscore", mode: params.publish_dir_mode

    input:
    file(bed_gz) // collected output from get_adj_zscore

    output:
    tuple file("adj_scores.bed.gz"), file("adj_scores.bed.gz.tbi")

    script:
    """
    merge_samples_adj_scores.py -r ./ | bgzip -c > adj_scores.bed.gz
    tabix -p bed adj_scores.bed.gz
    """
}


// Allele Balance Calculations

process collect_allele_counts {
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/${idSample}/collect_allele_counts", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(bam), file(bai)
    file(probes)
    file(ref)    
    file(fai)
    file(dict)

    output:
    tuple idPatient, idSample, file("${idSample}.allele_count.tsv")

    script:
    """
    gatk CollectAllelicCounts \
     -I $bam \
     -R $ref \
     -L $probes \
     -O ${idSample}.allele_count.tsv
    """
}

process agg_allele_counts{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/${idSample}/agg_allele_counts", mode: params.publish_dir_mode

    input:
    tuple idPatient, idSample, file(collected_allele_counts) // output from collect_allele_counts
    file(probes)

    output:
    tuple idPatient, idSample, file("${idSample}.agg.allele_count.bed")

    script:
    """
    cat $collected_allele_counts | awk '{print \$1"\t"\$2"\t"\$2+1"\t"\$3"\t"\$4"\t"\$5"\t"\$6}' > ${idSample}.allele_count.bed
    bgzip ${idSample}.allele_count.bed
    tabix -p bed ${idSample}.allele_count.bed.gz
    agg_single_allele_counts.py $probes ${idSample}.allele_count.bed.gz ${idSample} > ${idSample}.agg.allele_count.bed    
    """
}

process merge_all_allele_counts {
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/merge_all_allele_counts", mode: params.publish_dir_mode

    input:
    file(agged_allele_counts) // output from agg_allele_counts
    file(probes)

    output:
    tuple file("ab.sorted.tsv.gz"), file("ab.sorted.tsv.gz.tbi"), file("ab.header.tsv")
    
    script:
    """
    cat *.agg.allele_count.bed > all_samples.aggregate_probe_allele_counts.txt
    create_all_sample_allele_count_bed.py all_samples.aggregate_probe_allele_counts.txt $probes> ab.bed
    head -1 ab.bed > ab.header.tsv
    tail -n +2 ab.bed > ab.tsv
    bedtools sort -i ab.tsv > ab.sorted.tsv
    bgzip ab.sorted.tsv
    tabix -p bed ab.sorted.tsv.gz 
    """
}

process label_exons{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/label_exons", mode: params.publish_dir_mode

    input:
    file(exons_gz)

    output:
    tuple file("labeled_exons.bed.gz"), file("labeled_exons.bed.gz.tbi")

    script:
    """
    zcat $exons_gz | label_exon_number.py > labeled_exons.bed
    bgzip labeled_exons.bed
    tabix -p bed labeled_exons.bed.gz
    """
}

process combine_savvy_calls{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/Savvy_call_params", mode: params.publish_dir_mode

    input:
        file(beds)

    output:
        file("savvy_call_plotting_params.txt")

    script:
    """
    for i in *.bed; do
        cat \$i | awk -v VAR=\$i '{ gsub(".coverageBinner","",\$10); print \$1":"\$2"-"\$3"\t"\$4"\t"\$10"\t"VAR"\t"\$1"."\$2"-"\$3"\t"\$4}' >> tmp_savvy_call_plotting_params.txt
    done

    for i in {1..23}; do
        grep "^\$i" tmp_savvy_call_plotting_params.txt >> savvy_call_plotting_params.txt || echo "Not found"
    done

    touch thing.txt
    """
}

process testy2{
    label 'container_llab'
    label 'cpus_1'

    publishDir "${params.outdir}/CNV_Plotting/TESTAREA2", mode: params.publish_dir_mode

    input:
        tuple thing1, thing2, thing3, thing4

    output:
        file("testy.txt")

    """
    touch testy.txt
    echo $thing1 >> testy.txt
    echo "\n_\n"
    echo $thing2 >> testy.txt
    echo "\n_\n"
    echo $thing3 >> testy.txt
    echo "\n_\n"
    echo $thing4 >> testy.txt
    """

}

process cnv_plotter {
    label 'container_py3_pandas'
    label 'cpus_16'

    publishDir "${params.outdir}/CNV_Plotting/${idSample}/Plots", mode: params.publish_dir_mode

    input:
    tuple file(ab_bed_gz), file(ab_bed_gz_tbi), file(ab_header) // output from merge_all_allele_counts
    tuple file(adj_scores_bed_gz), file(adj_scores_bed_gz_tbi)
    // tuple all_caller_type,  all_vcf_idPatient, all_vcf_idSample, file(all_vcf), file(all_vcf_tbi) // should contain all vcfs in the named like ${idSample}.vcf
    file(all_vcf)
    tuple file(labeled_exons_bed_gz), file(labeled_exons_bed_gz_tbi) // output from label_exons
    tuple file(probe_cover_mean_stdev_bed_gz), file(probe_cover_mean_stdev_bed_gz_tbi), file(probe_cover_mean_stdev_bed)
    tuple region, sv_type, idSample, filename, region_without_colon
    file(beds)

    output:
    file("${idSample}.${region_without_colon}.${sv_type}.png")

    script:
    """
    for i in *.bed; do
        if [ "\$i" != "probe.cover.mean.stdev.bed" ]
        then 
            bedtools sort -i \${i} > \${i}.sorted 2>> error.txt
            bgzip \${i}.sorted 2>> error.txt
            tabix -p bed \${i}.sorted.gz 2>> error.txt
            echo "\${i}.sorted.gz" >> multiple_savvy_calls.txt 2>> error.txt
        fi
    done   


    cnv_plotter.py --sample $idSample \
    --vcf ${idSample}.vcf.gz \
    -o ${idSample}.${region_without_colon}.${sv_type}.png \
    --scores $adj_scores_bed_gz \
    --exons $labeled_exons_bed_gz \
    --window 100000 \
    --height 7 \
    --width 5 \
    --region $region \
    --title "${idSample} ${region} ${sv_type}" \
    --alt_allele_counts $ab_bed_gz \
    --alt_allele_headers $ab_header \
    --label_exons \
    --depth $probe_cover_mean_stdev_bed \
    --all_calls multiple_savvy_calls.txt 2>> error.txt

    """

}
















