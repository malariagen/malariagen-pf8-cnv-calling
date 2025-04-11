# Invoked by `run.sh`

PATH_TO_POSTPROCESSING: str = "../assets_pf8/03_execution/postprocessing"

# ===== imports =====
import pandas            as pd
import numpy             as np
import matplotlib.pyplot as plt

from matplotlib import gridspec

import os
import heapq
import sys


def _get_bounds(
    CALL_ID: str,
    CHROMOSOME: str,
    zoom: int = 30,
    breakpoint_buffer: int = 5000 # How much extra space you want to show beyond the breakpoint location
):
    """
    Helper function to help with deciding the x axis limits. Typically gets the length of the 
    call region then zooms out 30x so that the plot displays enough data either side of
    the call region to make accurate calls. 
    """
    gene_start = goi.loc[goi.CONTIG == CHROMOSOME, "START"].values[0]
    gene_end   = goi.loc[goi.CONTIG == CHROMOSOME, "END"].values[0]
    diff       = int(gene_end - gene_start)
    
    chromosome_length = chrom_lengths_dict[CHROMOSOME]
    
    # From Table 7 from https://www.malariagen.net/wp-content/uploads/2023/11/Pf7-supplementary-display-items-20221204.pdf
    if CHROMOSOME == "Pf3D7_08_v3": # HRP2 breakpoints can reach as low as 1373732
        lower_bound = np.min([1373732 - breakpoint_buffer, gene_start - diff * zoom])
        upper_bound = gene_end + diff * zoom
    elif CHROMOSOME == "Pf3D7_13_v3": # HRP3 breakpoints can reach as low as 2800004
        lower_bound = np.min([2800004 - breakpoint_buffer, gene_start - diff * zoom])
        upper_bound = gene_end + diff * zoom
    else:
        lower_bound = gene_start - diff * zoom
        upper_bound = gene_end + diff * zoom
    
    # ===== Ensuring plots don't escape the length of the corresponding chromosome =====
    lower_bound = np.max([0, lower_bound])
    upper_bound = np.min([chromosome_length, upper_bound])
    
    return gene_start, gene_end, lower_bound, upper_bound

def _assign_y_coords_to_df(df):
    """
    Helper function to assign a y value to horizontal tracks corresponding to genes
    from the GFF. This function ensures genes are evenly distributed without overlapping 
    whilst minimising the total vertical space consumed by subplot. 
    """
    genes = [(row.START, row.END, index) for index, row in df.iterrows()]
    genes.sort(key=lambda x: (x[0], x[1]))
    
    heap = []
    y_coords = [-1] * len(genes)
    
    for start, end, index in genes:
        if heap and heap[0][0] <= start:
            earliest_end, y = heapq.heappop(heap)
        else:
            y = len(heap)
        
        y_coords[index] = y
        heapq.heappush(heap, (end, y))
    
    df["Y_COORD"] = y_coords
    return df


def plot(
    SAMPLE_ID: str,
    CALL_ID  : str
):
    # ===== Checking input params =====
    assert SAMPLE_ID in os.listdir(PATH_TO_POSTPROCESSING), "Sample not found."
    assert CALL_ID in goi.CALL_ID.unique(), f"`CALL_ID` not found in `../assets_pf8/04_call_regions.tsv`. Choose from one of {goi.CALL_ID.unique()}"
    
    CHROMOSOME = goi.loc[goi.CALL_ID == CALL_ID, "CONTIG"].values[0]
    
    # ===== Loading pipeline outputs and wrangling =====
    sample_read_count_file = pd.read_csv(f"../assets_pf8/01_execution/readcounts/{SAMPLE_ID}.counts.tsv", sep = "\t", comment = "@")
    sample_read_count_file["INTERVAL_MIDPOINT"] = sample_read_count_file[["START", "END"]].mean(axis = 1)
    sample_read_count_file = sample_read_count_file.groupby("CONTIG")[["INTERVAL_MIDPOINT", "COUNT"]].apply(lambda s: ({
        "INTERVAL_MIDPOINT": np.array(s.INTERVAL_MIDPOINT),
        "COUNT": np.array(s.COUNT)
    })).to_dict()
    
    try:
        interval = pd.read_csv(os.path.join(PATH_TO_POSTPROCESSING, SAMPLE_ID, f"{SAMPLE_ID}.intervals.vcf.gz"), sep = "\t", skiprows = 24)
        interval.INFO = interval.INFO.apply(lambda x: x.replace("END=", "")).astype(int)
        interval["CN"] = interval.loc[:, interval.columns[-1]].apply(lambda x: x.split(":")[1]).astype(int)
        interval_dict = interval.groupby("#CHROM")[["POS", "INFO", "CN"]].apply(lambda s: [
            {"START": start, "END": end, "CN": cn}
            for start, end, cn
            in zip(s.POS.values, s.INFO.values, s.CN.values)
        ]).to_dict()
    except EOFError:
        print(f"EOFError for {SAMPLE_ID}, {CALL_ID}, interval.vcf.gz")
    
    try:
        segment = pd.read_csv(os.path.join(PATH_TO_POSTPROCESSING, SAMPLE_ID, f"{SAMPLE_ID}.segments.vcf.gz"), sep = "\t", skiprows = 35)
        segment.INFO = segment.INFO.apply(lambda x: x.replace("END=", "")).astype(int)
        segment["CN"] = segment.loc[:, segment.columns[-1]].apply(lambda x: x.split(":")[1]).astype(int)
        segment[["QA", "QS", "QSE", "QSS"]] = pd.DataFrame(
            segment[SAMPLE_ID].apply(lambda s: [int(v) for v in s.split(":")[3:]]).tolist(), 
            index = segment.index
        )
        segment_dict = segment.groupby("#CHROM").apply(lambda s: [
            {"START": start, "END": end, "CN": cn, "QA": qa, "QS": qs, "QSE": qse, "QSS": qss}
            for start, end, cn, qa, qs, qse, qss
            in zip(s.POS.values, s.INFO.values, s.CN.values, s.QA.values, s.QS.values, s.QSE.values, s.QSS.values)
        ]).to_dict()
    except EOFError:
        print(f"EOFError for {SAMPLE_ID}, {CALL_ID}, segments.vcf.gz")

    dCR = pd.read_csv(os.path.join(PATH_TO_POSTPROCESSING, SAMPLE_ID, f"{SAMPLE_ID}.dCR.tsv"), sep = "\t", skiprows = 18)
    dCR["INTERVAL_MIDPOINT"] = dCR[["START", "END"]].mean(axis = 1)

    for i, row in segment.iterrows():
        dCR.loc[
            (dCR.CONTIG == row["#CHROM"]) &
            (dCR.START >= row.POS) &
            (dCR.END <= row.INFO),
            "SEGMENT_CN"
        ] = row.CN
    
    dCR["DCR_SEGMENT_RESIDUAL"] = (dCR.LINEAR_COPY_RATIO - dCR.SEGMENT_CN).abs()
    
    assert not pd.isnull(dCR.DCR_SEGMENT_RESIDUAL).any(), "Not all residuals could be calculated. Ensure every dCR value is contained within one segment. "
    
    dCR_dict = dCR.groupby("CONTIG")[["INTERVAL_MIDPOINT", "LINEAR_COPY_RATIO", "DCR_SEGMENT_RESIDUAL"]].apply(lambda s: {
        "INTERVAL_MIDPOINT"   : np.array(s.INTERVAL_MIDPOINT),
        "LINEAR_COPY_RATIO"   : np.array(s.LINEAR_COPY_RATIO),
        "DCR_SEGMENT_RESIDUAL": np.array(s.DCR_SEGMENT_RESIDUAL),
    }).to_dict()
    
    # ===== Subsetting data objects by chromosome =====
    sample_read_count_file = sample_read_count_file[CHROMOSOME]
    dCR_dict               =               dCR_dict[CHROMOSOME]
    interval_dict          =          interval_dict[CHROMOSOME]
    segment_dict           =           segment_dict[CHROMOSOME]
    
    # ===== Initialising axes of figure =====
    fig = plt.figure(figsize = (7, 9))
    fig.suptitle(f"{SAMPLE_ID} - {CALL_ID}")
    gs = gridspec.GridSpec(11, 1, height_ratios = [1, 3, 3, 1, 1, 1, 1, 1, 0.5, 1.5, 1.5])
    
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1], sharex = ax0)
    ax2 = plt.subplot(gs[2], sharex = ax0)
    ax3 = plt.subplot(gs[3], sharex = ax0)
    ax4 = plt.subplot(gs[4], sharex = ax0)
    ax5 = plt.subplot(gs[5], sharex = ax0)
    ax6 = plt.subplot(gs[6], sharex = ax0)
    ax7 = plt.subplot(gs[7], sharex = ax0)
    
    ax8 = plt.subplot(gs[9])
    ax9 = plt.subplot(gs[10], sharex = ax8)
    
    plt.subplots_adjust(hspace = 0)

    # ===== Plotting gene annotations =====
    gene_start, gene_end, lower_bound, upper_bound = _get_bounds(CALL_ID, CHROMOSOME)

    gff_subset = gff.loc[(gff.CONTIG == CHROMOSOME) & (gff.END >= lower_bound) & (gff.START <= upper_bound)].reset_index(drop = True)
    gff_subset = _assign_y_coords_to_df(gff_subset)
    
    for index, row in gff_subset.iterrows():
        ax0.barh(
            y = row.Y_COORD, width = row.END - row.START,
            left = row.START, height = 0.9,
            align = "center", color = "black", alpha = 0.1)
        ax0.text(
            x = row.INTERVAL_MIDPOINT, y = row.Y_COORD, s = row.NAME,
            fontsize = 6, rotation = 45, va = "center", ha = "center")
    
    # ===== Masking and subsetting data =====
    dCR_dict_full = dCR_dict.copy() # Copying a fresh copy for later
    
    dCR_zoom_mask = (lower_bound <= dCR_dict["INTERVAL_MIDPOINT"]) & (dCR_dict["INTERVAL_MIDPOINT"] <= upper_bound)
    dCR_dict["INTERVAL_MIDPOINT"]    = dCR_dict["INTERVAL_MIDPOINT"][dCR_zoom_mask]
    dCR_dict["LINEAR_COPY_RATIO"]    = dCR_dict["LINEAR_COPY_RATIO"][dCR_zoom_mask]
    dCR_dict["DCR_SEGMENT_RESIDUAL"] = dCR_dict["DCR_SEGMENT_RESIDUAL"][dCR_zoom_mask]
    
    sample_read_count_file_full = sample_read_count_file.copy()
    readcounts_zoom_mask = (lower_bound <= sample_read_count_file["INTERVAL_MIDPOINT"]) & (sample_read_count_file["INTERVAL_MIDPOINT"] <= upper_bound)
    sample_read_count_file["INTERVAL_MIDPOINT"] = sample_read_count_file["INTERVAL_MIDPOINT"][readcounts_zoom_mask]
    sample_read_count_file["COUNT"]             = sample_read_count_file["COUNT"][readcounts_zoom_mask]
    
    # ===== Plotting CN data =====
    ax1.scatter(
        x = dCR_dict["INTERVAL_MIDPOINT"],
        y = dCR_dict["LINEAR_COPY_RATIO"],
        c = "r", s = 1, alpha = 0.7, zorder = 1)
    
    ax3.scatter(
        x = dCR_dict["INTERVAL_MIDPOINT"],
        y = dCR_dict["DCR_SEGMENT_RESIDUAL"],
        s = 0.6, alpha = 0.7, zorder = 1)
    
    ax8.scatter(
        x = dCR_dict_full["INTERVAL_MIDPOINT"],
        y = dCR_dict_full["LINEAR_COPY_RATIO"],
        c = "r", s = 0.2, alpha = 0.7, zorder = 1)
    
    for intvl in interval_dict:
        ax1.barh(
            y = intvl["CN"], width = intvl["END"] - intvl["START"],
            left = intvl["START"], height = 0.3, align = "center",
            color = "orange", alpha = 0.8, zorder = -1)
        ax8.barh(
            y = intvl["CN"], width = intvl["END"] - intvl["START"],
            left = intvl["START"], height = 0.4, align = "center",
            color = "orange", alpha = 0.8, zorder = -1)
    
    for seg in segment_dict:
        ax1.barh(
            y = seg["CN"], width = seg["END"] - seg["START"],
            left = seg["START"], height = 0.6, align = "center",
            color = "green", alpha = 0.3, zorder = -2)
        ax8.barh(
            y = seg["CN"], width = seg["END"] - seg["START"],
            left = seg["START"], height = 0.7, align = "center",
            color = "green", alpha = 0.3, zorder = -2)
    
    # ===== Plotting read counts =====
    ax2.bar(
        x      = sample_read_count_file["INTERVAL_MIDPOINT"],
        height = sample_read_count_file["COUNT"],
        width  = 500)
    
    ax9.bar(
        x      = sample_read_count_file_full["INTERVAL_MIDPOINT"],
        height = sample_read_count_file_full["COUNT"],
        width  = 500)

    # ===== Plotting QC metrics from segments.vcf.gz =====

    for metric, ax in zip(["QA", "QS", "QSS", "QSE"], [ax4, ax5, ax6, ax7]):
        # ----- Lots of gymnastics for handling placement of objects on figure carefully -----
        all_metric_values = [seg[ metric ] for seg in segment_dict]
        height = (max(all_metric_values) - min(all_metric_values)) / 5
        normalised_height = max([2, height])

        ax.set_yticks([])
        ax.set_ylabel(metric)

        if max(all_metric_values) - min(all_metric_values) < 10:
            ax.set_ylim(np.mean(all_metric_values) - 5, np.mean(all_metric_values) + 5)

        # From https://www.biorxiv.org/content/10.1101/2022.08.25.504851v1.full.pdf
        
        # QA measures the probability that all spanning intervals within the segment 
        # are consistent with the copy number assigned to the segment; 
        # QS represents the probability that at least one interval is consistent;
        # QSS and QSE represent the probability that the first and last intervals are consistent 
        # with being a change-point into and out of the segment copy number state, respectively. 

        # QS >= min(400,  max(50,   4 x NInt)) (for duplications)
        # QS >= min(1000, max(100, 10 x NInt)) (for CN=1 deletions)
        # QS >= min(1000, max(400, 10 x NInt)) (for CN=0 deletions)

        # Also this is interesting: https://gatk.broadinstitute.org/hc/en-us/articles/360035889891--Notebook-Correlate-gCNV-callset-metrics-and-annotations#2.3
        
        for seg in segment_dict:
            ax.barh(
                y = seg[metric], width = seg["END"] - seg["START"],
                left = seg["START"], height = normalised_height, align = "center",
                color = "gray", alpha = 0.3, zorder = -2)

            seg_left_edge_in  = (seg["START"] <= upper_bound) and (seg["START"] >= lower_bound)
            seg_right_edge_in = (seg["END"]   >= lower_bound) and (seg["END"]   <= upper_bound)
            seg_around_bounds = (seg["END"]   >= upper_bound) and (seg["START"] <= lower_bound)

            if metric == "QSS":
                ax.scatter(
                    x = seg["START"], y = seg[metric], marker = 9,
                    color = "green", zorder = -1, alpha = 0.2, s = 50)
            elif metric == "QSE":
                ax.scatter(
                    x = seg["END"], y = seg[metric], marker = 8,
                    color = "red", zorder = -1, alpha = 0.2, s = 50)
            else:
                ax.scatter(
                    x = seg["START"], y = seg[metric], marker = 9,
                    color = "blue", zorder = -1, alpha = 0.2, s = 50)
                ax.scatter(
                    x = seg["END"], y = seg[metric], marker = 8,
                    color = "blue", zorder = -1, alpha = 0.2, s = 50)
            
            if (seg_left_edge_in or seg_right_edge_in or seg_around_bounds):
                text_position = (seg["END"] + seg["START"]) / 2

                if (not (lower_bound <= text_position <= upper_bound)) and seg_around_bounds:
                    text_position = (lower_bound + upper_bound) / 2
                
                elif text_position >= upper_bound:
                    text_position = (upper_bound + seg["START"]) / 2
                
                elif text_position <= lower_bound:
                    text_position = (lower_bound + seg["END"]) / 2
                
                ax.text(
                    x = text_position, y = seg[metric],
                    s = f"{seg[metric]}.", fontweight = "bold",
                    ha = "center", va = "center", rotation = 15, zorder = 2)

    # ===== Fine-tuning plot settings =====
    # ----- Removing ticks from specific subplots -----
    plt.setp(ax0.get_xticklabels(), visible = False)
    plt.setp(ax1.get_xticklabels(), visible = False)
    plt.setp(ax2.get_xticklabels(), visible = False)
    plt.setp(ax3.get_xticklabels(), visible = False)
    plt.setp(ax4.get_xticklabels(), visible = False)
    plt.setp(ax5.get_xticklabels(), visible = False)
    plt.setp(ax6.get_xticklabels(), visible = False)
    ax7.tick_params(axis = "x", which = "both", labelbottom = True)
    plt.setp(ax8.get_xticklabels(), visible = False)
    ax9.tick_params(axis = "x", which = "both", labelbottom = True)
    
    ax0.set_yticks([])
    ax0.set_yticklabels([])

    # ----- Adding lines corresponding to CN ----
    ax1.grid(True, axis = "y")
    ax1.set_yticks([1, 2, 3, 4, 5])
    ax1.set_xlim(lower_bound, upper_bound)
    ax1.set_ylim(-0.05, 5.1)
    ax1.set_ylabel("COPY\nNUMBER")
    
    xlims = ax8.get_xlim()
    ax8.grid(True, axis = "y", linewidth = 0.3)
    ax8.set_yticks([1, 2, 3, 4, 5])
    ax8.set_xlim(xlims[0], xlims[1])
    ax8.set_ylim(-0.05, 5.1)
    ax8.set_ylabel("COPY\nNUMBER")
    
    ax8.axvspan(xlims[0], lower_bound, color = "black", alpha = 0.7)
    ax8.axvspan(upper_bound, xlims[1], color = "black", alpha = 0.7)
    
    # ----- Handling coverage axis lims -----
    y_lim = np.max((99, 2.5 * np.percentile(sample_read_count_file["COUNT"], 90))) # Upper y_lim is either 99 or 1.5 x 90th percentile
    if CALL_ID in ["HRP2", "HRP3"]:
        y_lim = 100
    ax2.set_ylim(0.1, y_lim)
    ax2.set_ylabel("COVERAGE")
    
    ax9.set_ylim(1, 499)
    ax9.set_ylabel("COVERAGE")
    
    ax9.axvspan(xlims[0], lower_bound, color = "black", alpha = 0.7)
    ax9.axvspan(upper_bound, xlims[1], color = "black", alpha = 0.7)
    
    # ----- Highlighting call region (e.g., gene of interest) ----
    for axis in [ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]:
        axis.axvline(x = gene_start, c = "red", alpha = 0.3)
        axis.axvline(x = gene_end,   c = "red", alpha = 0.3)

    # ----- Adding lines for different levels of residuals ----
    ax3.set_ylim(-0.1, 1.5)
    ax3.set_yticks([0.5, 1])
    
    ax3.axhline(y = 0,   c = "black", linewidth = 0.3)
    ax3.axhline(y = 0.5, c = "black", linewidth = 0.3)
    ax3.axhline(y = 1,   c = "black", linewidth = 0.3)
    ax3.axhline(y = 1.5, c = "black", linewidth = 0.3)
    ax3.set_ylabel("ABSOLUTE\nRESIDUALS")
    
    # ===== Finalising global plot settings =====
    plt.tight_layout()
    plt.subplots_adjust(hspace = 0)
    
    output_dir = f"results/{CALL_ID}/{SAMPLE_ID}/plot.png"
    
    os.makedirs(os.path.dirname(output_dir), exist_ok = True)
    
    plt.savefig(output_dir)
    plt.close("all")
    
    return


ANNOTATIONS_OF_INTEREST = ("protein_coding_gene", "ncRNA_gene", "pseudogene", "tRNA", "rRNA", "ncRNA", "snoRNA", "snRNA")

goi = pd.read_csv("../assets_pf8/04_call_regions.tsv", sep = "\t")

gff = pd.read_csv("../assets_pf8/PlasmoDB-55_Pfalciparum3D7.gff", sep = "\t", comment = "#",
                  names = ["CONTIG", "SOURCE", "TYPE", "START", "END", "SCORE", "STRAND", "PHASE", "INFO"])

gff = gff.loc[gff.TYPE.isin(ANNOTATIONS_OF_INTEREST) & gff.CONTIG.isin(goi.CONTIG.unique())]
gff["NAME"] = gff.TYPE.replace({"protein_coding_gene": "...", "ncRNA_gene": "ncRNA"})

# If the INFO field contains a gene `Name` value, then extract it and replace the `NAME` value
gff.loc[gff.INFO.str.contains("Name"), "NAME"] = (
    gff.loc[gff.INFO.str.contains("Name"), "INFO"]
    .apply(lambda x: [substr.replace("Name=", "") for substr in x.split(";") if "Name" in substr][0])
)

gff = gff[["CONTIG", "START", "END", "NAME"]].sort_values(["CONTIG", "START"]).reset_index(drop = True)
gff["INTERVAL_MIDPOINT"] = gff[["START", "END"]].mean(axis = 1)

# =======
chrom_lengths_dict = pd.read_csv("../assets_pf8/chromosomes_lengths.tsv", sep = "\t").set_index("chromosome").to_dict()["length"]
interval_list      = pd.read_csv("../assets_pf8/01_intervals_of_interest.interval_list",
                                 sep = "\t", comment = "@", names = ["CHROM", "START", "END", "_", "CONTIG"])

# Since the output folders are named after sample name, we can just do `os.listdir` to get list of samples
pf8_samples = os.listdir(PATH_TO_POSTPROCESSING)

LSB_JOBINDEX = int(sys.argv[1])
CHUNK_SIZE   = int(sys.argv[2])

first_index = (LSB_JOBINDEX - 1) * CHUNK_SIZE
last_index  = LSB_JOBINDEX * CHUNK_SIZE

samples = pf8_samples[first_index : last_index]

for sample in samples:
    plot(sample, "MDR1")
    plot(sample, "CRT")
    plot(sample, "HRP2")
    plot(sample, "GCH1")
    plot(sample, "HRP3")
    plot(sample, "PM2_PM3")