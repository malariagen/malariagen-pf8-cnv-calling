#!/usr/bin/env python

from optparse import OptionParser
import errno

import os
import numpy as np
import pandas as pd
import pysam
import collections

__version__ = "0.0.1"

if __name__ == '__main__':
    
    usage = 'usage: %prog [options]'
    description = "Run breakpoint read counts on a bam file"
    epilog = """
This runs the faceaway and same direction read counting around breakpoints pipeline on a bam file and creates a formatted output with sample name in first column, followed by various columns related to counts of different types of reads at different breakpoints

Examples:
    
    create_breakpoint_read_counts_20250128.py \
    --sample PH0906-C
    --bam_fn <INSERT PATH HERE>/74/PH0906-C/PH0906-C.bam \
    --tandem_dup_breakpoints_fn <INSERT PATH HERE>/66_faceaway_calls/20250128/tandem_dup_breakpoints_20250128.txt \
    --dup_trpinv_dup_breakpoints_fn <INSERT PATH HERE>/66_faceaway_calls/20250128/dup_trpinv_dup_breakpoints_20250128.txt \
    --region_size 600 \
    --region_offset 100 \
    --output_dir <INSERT PATH HERE>/66_faceaway_calls/20250128/samples
    
Version: {version}
""".format(version=__version__)

    OptionParser.format_epilog = lambda self, formatter: self.epilog
    parser = OptionParser(usage=usage, description=description, epilog=epilog)
    parser.add_option('-s', '--sample', dest='sample', help='sample ID [PH0906-C]', default='PH0906-C')
    parser.add_option('-i', '--bam_fn', dest='bam_fn', help='input bam filename [<INSERT PATH HERE>/74/PH0906-C/PH0906-C.bam]')
    parser.add_option('-t', '--tandem_dup_breakpoints_fn', dest='tandem_dup_breakpoints_fn', help='tandem duplication breakpoints filename [<INSERT PATH HERE>/66_faceaway_calls/20250128/tandem_dup_breakpoints_20250128.txt]')
    parser.add_option('-d', '--dup_trpinv_dup_breakpoints_fn', dest='dup_trpinv_dup_breakpoints_fn', help='DTD breakpoints filename [<INSERT PATH HERE>/66_faceaway_calls/20250128/dup_trpinv_dup_breakpoints_20250128.txt]')
    parser.add_option('-r', '--region_size', dest='region_size', type='int', help='Region size to search in base pairs [600]', default=600)
    parser.add_option('-f', '--region_offset', dest='region_offset', type='int', help='Region offset with respect to breakpoint in base pairs [100]', default=100)
    parser.add_option('-o', '--output_dir', dest='output_dir', help='output directory [<INSERT PATH HERE>/66_faceaway_calls/20250128/samples]')
    parser.add_option('-w', '--overwrite', action='store_true', dest='overwrite', default=False)
    options, args = parser.parse_args()
    
    try:
        
        def calc_proportion_faceaway(
            bam_fn,
            chrom,
            start_region, # 2-element list containing start and end positions (1-based) of region to right of first breakpoint in which to search for faceaway reads
            end_region,   # 2-element list containing start and end positions (1-based) of region to left of second breakpoint in which to search for faceaway reads
        ):
            samfile = pysam.AlignmentFile(bam_fn, "rb")
            num_faceaway = 0
            num_reads = 0
            iter = samfile.fetch(chrom, start_region[0], start_region[1])
            
            for x in iter:
                num_reads += 1
                if (
                    x.is_paired and
                    x.is_reverse and
                    (not x.mate_is_reverse) and
                    x.mpos >= end_region[0] and
                    x.mpos <= end_region[1] and
                    x.mapping_quality > 0 and # Important, as don't want ambiguously mapped reads
                    x.pos < x.mpos            # Important, e.g. in PfGCH1_promoter_dup_1, where start_region very close to end region
                ):
                    num_faceaway += 1
            
            if num_reads > 0:
                proportion_faceaway = num_faceaway/num_reads
            else:
                proportion_faceaway = float(np.nan)
            
            return(proportion_faceaway, num_faceaway, num_reads)
        
        # This looks for pairs of reads mapping in the same direction (indicative of an inversion)
        def calc_proportion_same_direction(
            bam_fn,
            chrom,
            start_region, # 2-element list containing start and end positions (1-based) of region to left of first breakpoint in which to search for faceaway reads
            end_region,   # 2-element list containing start and end positions (1-based) of region to left of second breakpoint in which to search for faceaway reads
            direction     # can be "forward" or "reverse"
        ):
            samfile = pysam.AlignmentFile(bam_fn, "rb")
            num_same_direction = 0
            num_reads = 0
            iter = samfile.fetch(chrom, start_region[0], start_region[1])
            
            for x in iter:
                num_reads += 1
                if direction=='forward':
                    if (x.is_paired and (not x.is_reverse) and (not x.mate_is_reverse) and x.mpos >= end_region[0] and x.mpos <= end_region[1]):
                        num_same_direction += 1
                if direction=='reverse':
                    if (x.is_paired and x.is_reverse and x.mate_is_reverse and x.mpos >= end_region[0] and x.mpos <= end_region[1]):
                        num_same_direction += 1
            
            if num_reads > 0:
                proportion_same_direction = num_same_direction/num_reads
            else:
                proportion_same_direction = float(np.nan)
            
            return(proportion_same_direction, num_same_direction, num_reads)

        def genotype_duplications(
            bam_fn,
            tandem_dup_breakpoints_fn,
            dup_trpinv_dup_breakpoints_fn,
            region_size                   = 600,
            region_offset                 = 100,
        ):
            df_tandem_dup_breakpoints = pd.read_csv(tandem_dup_breakpoints_fn, sep='\t')
            df_dup_trpinv_dup_breakpoints = pd.read_csv(dup_trpinv_dup_breakpoints_fn, sep='\t')
            
            genotype_results = collections.OrderedDict()
            genotype_results['faceaway'] = collections.OrderedDict()
            genotype_results['dup_trpinv_dup'] = collections.OrderedDict()
            
            for _, rec in df_tandem_dup_breakpoints.iterrows():
                genotype_results['faceaway'][rec.iloc[7]] = calc_proportion_faceaway(
                    bam_fn,
                    chrom = rec.iloc[0],
                    start_region = [
                        rec.iloc[2] - region_offset,
                        rec.iloc[2] - region_offset + region_size
                    ],
                    end_region = [
                        rec.iloc[3] - region_size,
                        rec.iloc[3],
                    ]
                )
            
            for _, rec in df_dup_trpinv_dup_breakpoints.iterrows():
                genotype_results['dup_trpinv_dup']["%s first" % rec.iloc[13]] = calc_proportion_same_direction(
                    bam_fn,
                    chrom = rec.iloc[0],
                    start_region = [
                        rec.iloc[1] - region_offset,
                        rec.iloc[1] - region_offset + region_size
                    ],
                    end_region = [
                        rec.iloc[3] - region_offset,
                        rec.iloc[3] - region_offset + region_size,
                    ],
                    direction = 'reverse'
                )
                genotype_results['dup_trpinv_dup']["%s second" % rec.iloc[13]] = calc_proportion_same_direction(
                    bam_fn,
                    chrom = rec.iloc[0],
                    start_region = [
                        rec.iloc[6] - region_size + region_offset,
                        rec.iloc[6] + region_offset
                    ],
                    end_region = [
                        rec.iloc[8] - region_size + region_offset,
                        rec.iloc[8] + region_offset,
                    ],
                    direction = 'forward'
                )
            
            return(genotype_results)

        def create_breakpoint_read_counts_file(
            sample,
            bam_fn,
            output_dir,
            tandem_dup_breakpoints_fn,
            dup_trpinv_dup_breakpoints_fn,
            region_size                   = 600,
            region_offset                 = 100,
            overwrite                     = False,
        ):
            
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                
            output_fn = f"{output_dir}/{sample}.tsv"
            
            if os.path.exists(output_fn) and not overwrite:
                print(f"File {output_fn} already exists, exiting")
        
            else:
                fo = open(output_fn, 'w')
            
                duplication_genotypes = genotype_duplications(
                    bam_fn,
                    tandem_dup_breakpoints_fn,
                    dup_trpinv_dup_breakpoints_fn
                )
        
                print("%s" % sample, end="", file=fo)
                for dup in duplication_genotypes['faceaway']:
                    results = duplication_genotypes['faceaway'][dup]
                    print("\t%d\t%d" % (results[1], results[2]), end="", file=fo)
                    
                for dup in duplication_genotypes['dup_trpinv_dup']:
                    results = duplication_genotypes['dup_trpinv_dup'][dup]
                    print("\t%d\t%d" % (results[1], results[2]), end="", file=fo)
                    
                print("\n", end="", file=fo)
        
                fo.close()

            pass
        
                
        create_breakpoint_read_counts_file(
            sample                        = options.sample,
            bam_fn                        = options.bam_fn,
            tandem_dup_breakpoints_fn     = options.tandem_dup_breakpoints_fn,
            dup_trpinv_dup_breakpoints_fn = options.dup_trpinv_dup_breakpoints_fn,
            region_size                   = options.region_size,
            region_offset                 = options.region_offset,
            output_dir                    = options.output_dir,
            overwrite                     = options.overwrite,
        )
    
    except IOError as e:
        if e.errno == errno.EPIPE:
            pass # ignore broken pipe
        else:
            raise
