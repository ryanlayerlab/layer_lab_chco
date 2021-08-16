#!/opt/miniconda/envs/chco/bin/python
import argparse
import utils
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pysam import VariantFile
import pandas as pd
import numpy as np
import pysam
import statistics
from random import randrange
from matplotlib.lines import Line2D

#{{{ def get_args():
def get_args():
    parser = argparse.ArgumentParser()


    parser.add_argument('--savvy_calls',
                        dest='savvy_calls',
                        help='BGZIPED and TABIXED Savvy calls set')

    parser.add_argument('--all_calls',
                        dest='all_calls',
                        help='file with path to BGZIPED and TABIXED Savvy calls set on each line')

    parser.add_argument('--scores',
                        dest='scores',
                        help='BGZIPED and TABIXED merged and adjusted scores')

    parser.add_argument('--exons',
                        dest='exons',
                        help='BGZIPED and TABIXED bed file with location of all exons')

    parser.add_argument('--vcf',
                        dest='vcf',
                        help='SNP/INDEL VCF"')

    parser.add_argument('--alt_allele_counts',
                        dest='alt_allele_counts',
                        help='tsv file with a normalized alterate allele count for each sample at each probe')

    parser.add_argument('--alt_allele_headers',
                        dest='alt_allele_headers',
                        help='headers for the tsv file with a normalized alterate allele count for each sample at each probe')

#    parser.add_argument('--cnvkit_calls',
#                        dest='cnvkit_calls',
#                        help='BGZIPED and TABIXED CNVKit calls set')
#
    parser.add_argument('--sample',
                        dest='sample',
                        help='Sample name')

    parser.add_argument('--region',
                        dest='region',
                        help='Target region')

    parser.add_argument('--window',
                        dest='window',
                        type=int,
                        default=50000,
                        help='Window (default 50000)')

    parser.add_argument('-o',
                        dest='outfile',
                        required=True,
                        help='Output file name')

    parser.add_argument('--legend_loc',
                        dest='legend_loc',
                        default='best',
                        help='Legend loctaion')
 
    parser.add_argument('--width',
                        dest='width',
                        type=float,
                        default=5,
                        help='Plot width (default 5)')

    parser.add_argument('--height',
                        dest='height',
                        type=float,
                        default=5,
                        help='Plot height (default 5)')

    parser.add_argument('--tick_line_length',
                        dest='tick_line_length',
                        type=float,
                        default=2,
                        help='Tick line width')

    parser.add_argument('--tick_line_width',
                       dest='tick_line_width',
                       type=float,
                       default=0.5,
                       help='Tick line width')

    parser.add_argument('--axis_line_width',
                       dest='axis_line_width',
                       type=float,
                       default=0.5,
                       help='Axis line width')

    parser.add_argument('--axis_label_size',
                        dest='axis_label_size',
                        type=int,
                        default=8,
                        help='Axis label font size')

    parser.add_argument('--tick_label_size',
                        dest='tick_label_size',
                        type=int,
                        default=8,
                        help='Axis tick label font size')

    parser.add_argument('--x_label',
                        dest='x_label',
                        help='X axis label')

    parser.add_argument('--y_label',
                        dest='y_label',
                        help='Y axis label')

    parser.add_argument("--title",
                        dest="title",
                        help="Plot title (title or title;size;location)")

    parser.add_argument("--depth",
                        dest="depth",
                        help="file with depth/rpm data")

    parser.add_argument('--label_exons',dest='label_exons', action='store_true',help="set to label each exon with it's gene name and exon number. requires that the file in --exons has gene and exon number listed in the two fields immediatelty after teh normal bed format items")

    args = parser.parse_args()

    return args
#}}}


def get_file_name(f):
    """
    removes file endings and file path information
    example: '/my/file/path.txt' -> 'path'
    f: str, file to have
    """
    return f.split('/')[-1].split('.')[0]

def get_region_alt_counts(chr,st,fin,bed_tbx,sample_col):
    res = bed_tbx.fetch(chr,st,fin)
    res_dict = {'chr':[],'start':[],'end':[],'avg':[],'sample':[],'zscore':[],'std':[]}
    for r in res:
        row = r.strip().split('\t')
        vals = [float(x) for x in row[4:]]
        mean = sum( vals) / len(vals)
        std = statistics.stdev(vals)
        if std != 0:
            zscore = (float(row[sample_col]) - mean) / std
        else:
            zscore = None
        res_dict['chr'].append(row[0])
        res_dict['start'].append(row[1])
        res_dict['end'].append(row[2])
        res_dict['avg'].append(mean)
        res_dict['sample'].append(row[sample_col])
        res_dict['zscore'].append(zscore)
        res_dict['std'].append(std)
    return pd.DataFrame(res_dict)

def format_axis(ax, args):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['bottom'].set_linewidth(args.axis_line_width)
    ax.spines['left'].set_linewidth(args.axis_line_width)
    ax.spines['right'].set_linewidth(args.axis_line_width)
    ax.tick_params(axis='both',
                   which='major',
                   labelsize=args.axis_label_size,
                   width=args.tick_line_width,
                   length=args.tick_line_length)

def mark_intervals(ax, intervals, color,a=0.25,ymin=0,ymax=1):
    for interval in intervals:
        ax.axvspan(max(ax.get_xlim()[0], interval.start),
                   min(ax.get_xlim()[1], interval.end),
                   ymin,ymax,
                   alpha=a, color=color)

def label_exons(ax, intervals):
    texts = []
    for interval in intervals:
        texts.append(ax.text(interval.start, 1, interval.data[0] + ' ' + interval.data[1], fontsize=3,rotation=90))
    return texts

def get_savvy_calls(savvy_bed_file, target_sample, target_interval):
    calls = utils.get_intervals_in_region(target_interval, savvy_bed_file)
    sample_calls = []
    for call in calls:
        curr_sample = call.data[6].split('.')[0]
        if curr_sample == target_sample:
            sample_calls.append(call)
    return(sample_calls)

def get_all_savvy_calls_in_target(savvy_bed_file, target_interval, ignore=None):
    calls = utils.get_intervals_in_region(target_interval, savvy_bed_file,ignore=ignore)
    sample_calls = []
    for call in calls:
        curr_sample = call.data[6].split('.')[0]
        sample_calls.append(call)
    return(sample_calls)

def get_cnvkit_calls(cnvkit_bed_file, target_sample, target_interval):
    calls = utils.get_intervals_in_region(target_interval, cnvkit_bed_file)
    sample_calls = []
    for call in calls:
        cn = int(call.data[1]) 
        if cn == 2:
            continue

        sample = call.data[5]

        if sample != target_sample.split('_')[0]:
            continue

        sample_calls.append(call)
    return(sample_calls)

def man_plot(axs, infile, target):
    cmap= plt.get_cmap('tab10')
    sets = []
    i = 0
    chrm = None
    curr_set = [[],[],[]]
    chrm_labels = [[],[]]
    xs=[]
    means=[]
    stds=[]
    for l in open(infile):
        A = l.rstrip().split('\t')
        if A[0] != target.chrom: continue
        xs.append(int(A[1]))
        means.append(float(A[4]))
        stds.append(float(A[5]))

    axs[0].scatter(xs,means,s=.1,color='#5683d6')
    axs[1].scatter(xs,stds,s=.1,color='#5683d6')

    axs[0].set_ylabel('Mean', fontsize=6)

    axs[1].set_ylabel('Stdev', fontsize=6)

    axs[1].set_xticks(chrm_labels[0])
    axs[1].set_xticklabels(labels=chrm_labels[1])
    axs[1].axvspan(target.start, target.end,alpha=.5, color='#ff5025')
    axs[0].axvspan(target.start, target.end,alpha=.5, color='#ff5025')

    return axs


def main():

    args = get_args()
    call_colors = ['#37dc94', '#8931EF', '#F2CA19', '#FF00BD', '#0057E9', '#87E911', '#E11845']

    chrom=args.region.split(':')[0]

    target_region = utils.Interval(
            chrom=chrom,
            start=max(0,
                      int(args.region.split(':')[1].split('-')[0]) - args.window),
            end=int(args.region.split(':')[1].split('-')[1]) + args.window,
            data=None)


    samples = utils.get_header(args.scores)[0].split('\t')[4:]
    sample_i = samples.index(args.sample)
    scores = utils.get_intervals_in_region(target_region,
                                           args.scores)

    exons = utils.get_intervals_in_region(target_region,
                                           args.exons)

    regions_pos = []
    means = []
    stdevs = []
    xs = []
    obs_hets = []
    mean_hets = []
    obs_homs = []
    mean_homs = []
    other = []
    alt_xs = []
    alt_ys = []
    alt_stds = []
    alt_means = []

    non_samp_xs = []    
    non_samp_pos = []
    vcf = VariantFile(args.vcf)

    for exon in scores:
        scores = [float(x) for x in exon.data[1:]]
        regions_pos.append(exon.end)
        means.append(np.mean(scores))
        stdevs.append(np.std(scores))
        xs.append(scores[sample_i])
        non_samp_xs += scores
        non_samp_pos += [exon.end] * len(scores)
        gt_stat =  utils.get_gt_stats(exon, target_sample=args.sample, vcf_handle=vcf)
        if gt_stat:
            obs_homs.append((exon.end,gt_stat[0][1]))
            obs_hets.append((exon.end,gt_stat[1][1]))

            mean_homs.append((exon.end,gt_stat[0][2]))
            mean_hets.append((exon.end,gt_stat[1][2]))

    plt.rcParams.update({'font.size': 6})
    #fig, axs = plt.subplots(4,1,gridspec_kw={'height_ratios': [3, 1, 1, 1]})
    axs = []
    fig = plt.figure(constrained_layout=True)
    gs = fig.add_gridspec(20,1, hspace=1.0, wspace=1.0)
    gs1 = fig.add_gridspec(20, 1, hspace=0.05, wspace=0.05)
    axs.append(fig.add_subplot(gs[:6,:]))
    axs.append(fig.add_subplot(gs[6:12,:]))
    axs.append(fig.add_subplot(gs[12:14,:]))
    axs.append(fig.add_subplot(gs1[15:17,:]))
    axs.append(fig.add_subplot(gs1[17:19,:]))
    
    fig.tight_layout()
    for a in axs:
        format_axis(a,args)
    ax = axs[0]
    #fig = plt.figure(figsize=(args.width,args.height), dpi=300)
    #rows=1
    #cols=1
    #outer_grid = gridspec.GridSpec(rows,
    #                               cols,
    #                               wspace=0.0,
    #                               hspace=0.25)
    #inner_grid = gridspec.GridSpecFromSubplotSpec(\
    #        1,
    #        1,
    #        subplot_spec=outer_grid[0],
    #        wspace=0.0,
    #        hspace=0.0)
    #ax = fig.add_subplot(inner_grid[0])
    if args.title:
        if ';' in args.title:
            text, size, loc = args.title.split(';')
            fontdict = {'fontsize':int(size)}
            ax.set_title(text, fontdict=fontdict, loc=loc)
        else:
            ax.set_title(args.title, loc='left')

    if args.alt_allele_counts and args.alt_allele_headers:
        tbx = pysam.TabixFile(args.alt_allele_counts)
        for line in open(args.alt_allele_headers):
            t_headers = line.strip().split('\t')
            break
        samp = t_headers.index(args.sample)
        region = args.region.split(':')
        c = region[0]
        s = max(0, int(region[1].split('-')[0]) - args.window)
        e = int(region[1].split('-')[1]) + args.window
        alt_info =  get_region_alt_counts(c,s,e,tbx,samp)
        alt_xs = list(alt_info['start'])
        alt_ys = list(alt_info['zscore'])
        alt_means = list(alt_info['avg'])
        alt_stds = list(alt_info['std'])
        alt_score = list( float(x) for x in alt_info['sample'])

    lns2b = ax.plot(non_samp_pos,
                   non_samp_xs,
                   '-o',
                   lw = 0,
                   markersize=1,
                   label='Other Sample Coverage',
                   c='#e2e2e2', alpha=.1)

    lns1 = ax.plot(regions_pos,
                   means,
                   marker='_',
                   markersize=3,
                   lw = 0,
                   label='Mean Pop. Coverage',
                   c='#5683d6')

    for i in range(len(means)):
        ax.plot([regions_pos[i],regions_pos[i]],
                [means[i]-stdevs[i],means[i]+stdevs[i]],
                lw=1,
                color='#5683d6',
                alpha=0.5)


    lns2 = ax.plot(regions_pos,
                   xs,
                   '-o',
                   lw = 0,
                   markersize=1,
                   label='Sample Coverage',
                   c='#5683d6')
    #lns3 = ax.plot([het[0] for het in obs_hets],
    #               [het[1] for het in obs_hets],
    #               lw = 0,
    #               marker='+',
    #               markersize=2,
    #               mew=0.25,
    #               label='Obs het',
    #               mec='black')

    #lns4 = ax.plot([het[0] for het in mean_hets],
    #               [het[1] for het in mean_hets],
    #               lw = 0,
    #               marker='+',
    #               markersize=2,
    #               mew=0.25,
    #               label='Exp het',
    #               mec='green')


    #lns5 = ax.plot([hom[0] for hom in obs_homs],
    #               [hom[1] for hom in obs_homs],
    #               lw = 0,
    #               marker='x',
    #               markersize=2,
    #               label='Obs hom',
    #               mew=0.25,
    #               mec='black')

    #lns6 = ax.plot([hom[0] for hom in mean_homs],
    #               [hom[1] for hom in mean_homs],
    #               lw = 0,
    #               marker='x',
    #               markersize=2,
    #               label='Exp hom',
    #               mew=0.25,
    #               mec='green')

    ax.spines['left'].set_visible(False)

    #ax2 = ax.twinx()
    ax2 = axs[1]
    #lns8 = ax2.fill_between([ int(x) for x in alt_xs], [ alt_means[i] - alt_stds[i] for i in range(len(alt_ys))],
    #                       [ alt_means[i] + alt_stds[i] for i in range(len(alt_ys))],
    #                       alpha=0.2,color='green' )

    lns8 = ax2.plot([ int(x) for x in alt_xs],
                   alt_means,
                   marker='_',
                   markersize=3,
                   lw = 0,
                   label='Mean Pop. Alt. Count',
                   c='#ff5025')

    alt_xs = [int(x) for x in alt_xs]
    for i in range(len(alt_xs)):
        ax2.plot([alt_xs[i],alt_xs[i]],
                [alt_means[i]-alt_stds[i],alt_means[i]+alt_stds[i]],
                lw=1,
                color='#ff5025',
                alpha=0.5)

    #lns10 = ax2.scatter([ int(x) for x in alt_xs], alt_score, label='Norm Alt Count', s=1, c='black')

    lns10 = ax2.plot(alt_xs,
                   alt_score,
                   '-o',
                   lw = 0,
                   markersize=1,
                   label='Sample Pop. Alt. Count',
                   c='#ff5025')

    format_axis(ax, args)

    ax.set_ylabel('Z-score', fontsize=6)

    #ax.set_xlabel(str(args.region.split(':')[0]), fontsize=args.axis_label_size)

    xmin,xmax = [int(x) for x in ax.get_xlim()]
    target_interval = utils.Interval(chrom=chrom,
                                     start=max(0,xmin),
                                     end=xmax,
                                     data=None)

    if args.savvy_calls and args.sample:
        xmin,xmax = [int(x) for x in ax.get_xlim()]
        target_interval = utils.Interval(chrom=chrom,
                                         start=max(0,xmin),
                                         end=xmax,
                                         data=None)
        savvy_calls = get_savvy_calls(args.savvy_calls,
                                      args.sample,
                                      target_interval)
        all_savvy_calls = get_all_savvy_calls_in_target(args.savvy_calls,target_interval)
        mark_intervals(ax, savvy_calls, call_colors[0],a=0.2)
        mark_intervals(ax2, savvy_calls, call_colors[0],a=0.2)
    if args.all_calls:
        all_calls = []
        number_of_calls_to_plot = 0
        for line in open(args.all_calls,'r'):
            line = line.strip()
            print(line)
            try:
                tbx_calls = pysam.TabixFile(line)
                raw_calls = tbx_calls.fetch(target_interval.chrom, target_interval.start, target_interval.end)
                calls = []
                for l in raw_calls:
                    if args.sample not in l: continue
                    row = l.strip().split('\t')
                    calls.append(utils.Interval(chrom=row[0],
                                         start=int(row[1]),
                                         end=int(row[2]),
                                         data=None)) 
            except ValueError:
                print('ValueError with file ' + line)
                calls = []
            #if there is 1 or more calls, add 1 to the count
            if len(calls) > 0:
                #number_of_calls_to_plot += 1
                number_of_calls_to_plot += len(calls)
            all_calls.append(calls)
            print(calls)
        box_num = 0
        for i,calls in enumerate(all_calls):
            if len(calls) == 0: continue
            for call in calls:
                if number_of_calls_to_plot == 0: break
                ymin = 1.0 / number_of_calls_to_plot * box_num
                ymax = ymin + (1.0 / number_of_calls_to_plot)
                mark_intervals(ax, [call], call_colors[i],a=0.2,ymin=ymin,ymax=ymax)
                mark_intervals(ax2, [call], call_colors[i],a=0.2,ymin=ymin,ymax=ymax)
                box_num += 1

    if args.exons:
        probes = []
        for line in open(args.depth,'r'):
            row = line.strip().split('\t')
            if str(row[0]) != str(target_region.chrom):
                continue
            probes.append(utils.Interval(chrom=row[0],start=int(row[1]),end=int(row[2]),data=None))
        tick_xs = [x.start for x in exons]
        tick_labs = [' '.join(x.data) for x in exons]
        #ax2.set_xticks(tick_xs)
        #ax2.set_xticklabels(tick_labs,fontsize= 3,)
        axs[2].set_xticks(tick_xs)
        axs[2].set_xticklabels(tick_labs,fontsize= 3,)
        #density_of_calls_at_each_probe = [ len(get_all_savvy_calls_in_target(args.savvy_calls,x)) for x in probes]
        #x_vals = [ x.start for x in probes]
        #indexes = [i for i,x in enumerate(density_of_calls_at_each_probe) if x > 0]
        #x_vals = [ x_vals[i] for i in indexes]
        #y_vals = [ density_of_calls_at_each_probe[i] for i in indexes]
        #axs[2].scatter(x_vals, y_vals, s=.1, color=call_colors[0], marker='*')
        legend_elements = []
        #legend_elements.append(Line2D([0], [0],
        #              marker='o',
        #              color=call_colors[0],
        #              label=get_file_name(args.savvy_calls),
        #              markerfacecolor=call_colors[0],
        #              markersize=2,
        #              linewidth=0))
        colors = call_colors
        markers = ['o', '^', 's', 'P', 'D', '*']
        for i,line in enumerate(open(args.all_calls,'r')):
            f = line.strip()
            try:
                density_of_calls_at_each_probe = [ len(get_all_savvy_calls_in_target(f,x,ignore=args.sample)) for x in probes]
            except ValueError:
                print('ValueError in file ' + f)
                density_of_calls_at_each_probe = [0 for x in probes]
            x_vals = [ x.start for x in probes]
            indexes = [j for j,x in enumerate(density_of_calls_at_each_probe) if x > 0]
            x_vals = [ x_vals[j] for j in indexes]
            y_vals = [ density_of_calls_at_each_probe[i] for i in indexes]
            axs[2].scatter(x_vals, y_vals,s=.1,color=colors[i],marker='*')
            legend_elements.append(Line2D([0], [0],
                              marker='o',
                              color=colors[i],
                              label=get_file_name(f),
                              markerfacecolor=colors[i],
                              markersize=2,
                              linewidth=0))
        axs[1].set_xticks([])
        axs[1].set_xticklabels([])
        #axs[2].set_xticks([])
        #axs[2].set_xticklabels([])
        axs[3].set_xticks([])
        axs[3].set_xticklabels([])

        axs[2].set_ylabel('Num. Calls',fontsize=6)
        axs[2].tick_params(bottom=False)
        #axs[2].axvspan(target_region.start, target_region.end,alpha=.5, color='#ff5025')
        axs[2].legend(handles=legend_elements,frameon=False,prop={'size': 5},bbox_to_anchor=(1.10,1), borderaxespad=0)
        if max(density_of_calls_at_each_probe) < 3 and False:
            print('Changed!!!')
            ticks = list(range(max(density_of_calls_at_each_probe)+1))
            print(ticks)
            axs[2].set_yticks(ticks)
            axs[2].set_yticklabels(ticks)
        else:
            print('no change!!')
            print(axs[2].get_yticks())
            print(axs[2].get_yticklabels())
#    if args.cnvkit_calls and args.sample_name:
#        xmin,xmax = [int(x) for x in ax.get_xlim()]
#        target_interval = utils.Interval(chrom=chrom,
#                                         start=xmin,
#                                         end=xmax,
#                                         data=None)
#        cnvkit_calls = get_cnvkit_calls(args.cnvkit_calls,
#                                        args.sample_name,
#                                        target_interval)
#
#        mark_intervals(ax, cnvkit_calls, 'blue')
#
#
    #lns = lns1 + lns2 + lns3 + lns4 + lns5 + lns6 + lns7
    lns = lns1 + lns2 + lns8 + lns10 + lns2b
    
    man_plot([axs[3],axs[4]], args.depth, target_region)
    ax2.set_ylabel('Normalized Alt. Allele Count', fontsize=6)
    labs = [l.get_label() for l in lns]
    ax.set_ylim([-7, 7])
    leg = ax.legend(lns,
                    labs,
                    frameon=False,
                    prop={'size': 5},bbox_to_anchor=(1.10,1), borderaxespad=0)
    for i in  range(0,len(axs)):
        axs[i].tick_params(axis='both', which='major', labelsize=4,length=2,width=.5)
        axs[i].tick_params(axis='both', which='major', labelsize=4,length=2,width=.5)

    #ax.tick_params(axis='both', which='major', labelsize=4,length=2,width=.5)
    #ax.tick_params(axis='both', which='minor', labelsize=4,length=2,width=.5)
    ax2.tick_params(axis='both', which='major', labelsize=4,length=2,width=.5)
    ax2.tick_params(axis='both', which='minor', labelsize=4,length=2,width=.5)
    axs[1].tick_params(axis='both', which='major', labelsize=4,length=2,width=.5)
    axs[1].tick_params(axis='both', which='minor', labelsize=4,length=2,width=.5)

    axs[2].tick_params(axis='both', which='major', labelsize=4,length=2,width=.5)
    axs[2].tick_params(axis='both', which='minor', labelsize=4,length=2,width=.5)

    axs[3].set_xticks([])

    #ax2.tick_params(axis='x', which='major', labelsize=2, rotation=90)
    #ax2.tick_params(axis='x', which='minor', labelsize=2, rotation=90)
    axs[2].tick_params(axis='x', which='major', labelsize=2, rotation=90)
    axs[2].tick_params(axis='x', which='minor', labelsize=2, rotation=90)

    for i in [0,1,2,3]:
        axs[i].spines['top'].set_visible(False)
        axs[i].spines['right'].set_visible(False)
        axs[i].spines['bottom'].set_visible(False)
        axs[i].spines['left'].set_visible(False)
    axs[1].spines['left'].set_visible(False)
    axs[-1].spines['left'].set_visible(False)
    axs[-1].set_xlabel('Chromosome ' + str(target_region.chrom), fontsize=6)
    axs[2].spines['bottom'].set_visible(True)

    ax.set_xticks([])
    ax.set_xticklabels([])

    regional_min = min(ax.get_xlim()[0],ax2.get_xlim()[0])
    regional_max = max(ax.get_xlim()[1],ax2.get_xlim()[1])
    
    ax.set_xlim((regional_min,regional_max))
    ax2.set_xlim((regional_min,regional_max))
    
    mins = []
    maxs = []
    
    for i in range(2,5):
        mins.append(axs[i].get_xlim()[0])
        maxs.append(axs[i].get_xlim()[1])
    for i in range(2,5):
        axs[i].set_xlim((min(mins),max(maxs)))

    axs[2].set_xlim((regional_min,regional_max))

    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.gca().spines['left'].set_color('none')
    plt.savefig(args.outfile,bbox_inches='tight',dpi=600)

if __name__ == '__main__': main()
