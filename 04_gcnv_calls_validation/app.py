BROAD_CN_RATIO_PATH  = "<INSERT PATH HERE>/data/cnv-broad-results/samples/"
DIAGNOSTIC_PLOT_PATH = "<INSERT PATH HERE>/malariagen-pf8-cnv-calling/04_gcnv_calls_validation/results"
PIPELINE_OUTPUT_PATH = "<INSERT PATH HERE>/malariagen-pf8-cnv-calling/assets_pf8/03_execution"
SAMPLE_META_PATH     = "<INSERT PATH HERE>/malariagen-pf8-cnv-calling/assets_pf8/Pf_8_samples_20241212.txt"

# ================================================================================

import json
import os

import streamlit            as st
import pandas               as pd
import matplotlib.pyplot    as plt
import numpy                as np
import plotly.graph_objects as go

# ================================================================================

@st.cache_data
def cache_load_plot_paths(CALL_ID, RAW_INPUT_STR):
    """
    Caches a map of {sample id : path to corresponding figure}. 
    
    Users can optionally provide Python lists containing specific samples of interest to display only those samples. 
    """

    BASE_PATHS_FOR_FIGURES = f"{DIAGNOSTIC_PLOT_PATH}/{CALL_ID}"
    
    if RAW_INPUT_STR == "":
        return {
            sample: os.path.join(BASE_PATHS_FOR_FIGURES, sample, "plot.png")
            for sample in os.listdir(BASE_PATHS_FOR_FIGURES)
        }
        
    else:
        l_samples = eval(RAW_INPUT_STR)
        return {
            sample: os.path.join(BASE_PATHS_FOR_FIGURES, sample, "plot.png")
            for sample in os.listdir(BASE_PATHS_FOR_FIGURES) 
            if sample in l_samples
        }




@st.cache_data
def cache_load_cnv_calls():
    """
    Load in the Broad Institute's old CNV calls as well as Pf7's drug resistance file for concordance analysis.
    """
    interim_1 = pd.read_csv(
        "app_files/Pf7_duplication_calls_20210114.tsv", sep = "\t"
    ).set_index("Sample")
    interim_1.columns = ["Interim1:"+col.replace("PM2", "PM2_PM3") for col in interim_1.columns.str.upper()]
    interim_1 = interim_1[[col for col in interim_1.columns if col.endswith("_CALL") or col.endswith("_ONLY")]]

    interim_2 = pd.read_csv(
        "app_files/Pf7_duplication_calls_20210416.tsv", sep = "\t"
    ).set_index("Sample")
    interim_2.columns = ["Interim2:"+col.replace("PM2", "PM2_PM3") for col in interim_2.columns.str.upper()]
    interim_2 = interim_2[[col for col in interim_2.columns if col.endswith("_CALL") or col.endswith("_ONLY")]]

    pf7_resistance = pd.read_csv(
        "app_files/Pf7_inferred_resistance_status_classification.tsv", sep = "\t"
    ).set_index("Sample")[
        ["Mefloquine", "Piperaquine", "HRP2", "HRP3"]
    ].applymap(lambda val: {
        "Sensitive"   :  0,
        "Resistant"   :  1,
        "Undetermined": -1,
        "nodel"       :  0,
        "del"         :  1,
        "uncallable"  : -1
        }[val]
    ).rename(columns = {
        "Mefloquine" : "Pf7: MDR1",
        "Piperaquine": "Pf7: PM2_PM3",
        "HRP2"       : "Pf7: HRP2",
        "HRP3"       : "Pf7: HRP3"
    })

    draft_calls = pd.read_csv("app_files/draft_coverage_calls.tsv", sep = "\t").set_index("SAMPLE")
    draft_calls.columns = "Pf8: " + draft_calls.columns

    pf7_qc_pass_samples_list = pf7_resistance.index.tolist()

    Pf7_resistance_calls = pd.concat([interim_1, interim_2, pf7_resistance, draft_calls], axis = 1).fillna("").map(lambda x: int(x) if x != "" else x)

    pf8_meta = pd.read_csv(SAMPLE_META_PATH, sep = "\t", usecols = ["Sample", "QC pass", "Sample was in Pf7"]).set_index("Sample")

    Pf8_calls = Pf7_resistance_calls.reset_index().rename(columns = {"index": "Sample"})
    Pf8_calls = Pf8_calls.loc[
        Pf8_calls.Sample.isin(pf8_meta.loc[pf8_meta["QC pass"] == True].index.values)
    ].reset_index(drop = True).replace({"": "no data"})

    return Pf7_resistance_calls, pf7_qc_pass_samples_list, Pf8_calls




@st.cache_data
def cache_load_bounds_assets():
    """
    Loading in call regions (e.g., CRT + CRT start and end coordinates) and information
    on length of each chromosome to assist with plotting. 
    """
    goi = pd.read_csv(
        "../assets_pf8/04_call_regions.tsv", sep = "\t")

    chrom_lengths = pd.read_csv(
        "../assets_pf8/chromosomes_lengths.tsv",
        sep = "\t").set_index("chromosome").to_dict()["length"]
    
    return goi, chrom_lengths




@st.cache_data
def cache_load_contig_ploidy_fails():
    """
    Some contigs will not have a ploidy of 1, so these are not suitable for making
    CNV calls with. Gathering all contigs where this is the case. 
    """
    failed_train_ploidy = pd.read_csv(
        f"{PIPELINE_OUTPUT_PATH}/failed_training_ploidy.tsv",
        sep = "\t", header = None, names = ["SAMPLE_ID", "CHROM", "PLOIDY", "GQ"])
    
    failed_test_ploidy = pd.read_csv(
        f"{PIPELINE_OUTPUT_PATH}/failed_test_ploidy.tsv",
        sep = "\t", header = None, names = ["SAMPLE_ID", "CHROM", "PLOIDY", "GQ"])
    
    return pd.concat([failed_train_ploidy, failed_test_ploidy]).reset_index(drop = True)




def update_sample_index():
    """
    Helper function for allowing the "Next" and "Previous" buttons to track sample index
    when clicking through plots. 
    """
    st.session_state.sample_index = SAMPLE_TO_INDEX_DICT[st.session_state.sample_selectbox]




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
    gene_start = goi.loc[goi.CALL_ID == CALL_ID, "START"].values[0]
    gene_end   = goi.loc[goi.CALL_ID == CALL_ID, "END"].values[0]
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




# ===== Set up page format =====
st.set_page_config(layout = "wide")

sidebar, main = st.columns([2, 5])

browser_mode = sidebar.radio(
    "Select browser mode",
    ["General browsing", "Modify draft calls", "Modification stats", "Inspect modifications", "Inspect final calls"])

sidebar.divider()

# ----- Quality of life features for sample filtering ----
Pf7_resistance_calls, pf7_qc_pass_samples_list, Pf8_calls = cache_load_cnv_calls()

if browser_mode == "General browsing":
    # ----- Initialise text input widget -----
    RAW_INPUT_STR = sidebar.text_input("""Provide some Python code that provide the sample names - e.g., `["FP0008-C", "FP0009-C"]`""")
    
    # ----- Allow user selection of `CALL ID` -----
    CALL_IDS = os.listdir(DIAGNOSTIC_PLOT_PATH)
    CALL_ID  = sidebar.selectbox("Call ID:", CALL_IDS)
    
    if sidebar.checkbox("16203 Pf7 samples only"):
        if CALL_ID in ["CRT"]:
            main.error(f"{CALL_ID} was not called in Pf7")
            st.stop()
        
        RAW_INPUT_STR = str(pf7_qc_pass_samples_list)
    
    sidebar.divider()




if browser_mode == "Modify draft calls":
    # Allow user to choose call region, then display confusion matrix dataframe
    # When the user selects a row, show all samples belonging to that genotype combination
    # e.g., all samples that were GT=1 in Pf7 but GT=0 in Pf8. 

    # Further down in this script, there is code that allows the user to modify the call. 
    n_total_samples = len(Pf8_calls)

    CALL_IDS = os.listdir(DIAGNOSTIC_PLOT_PATH)
    CALL_ID  = sidebar.selectbox("Call ID:", CALL_IDS)
    
    if CALL_ID in ["MDR1", "HRP2", "HRP3", "PM2_PM3"]:
        groupby = Pf8_calls.groupby([f"Pf7: {CALL_ID}", f"Pf8: {CALL_ID}"]).apply(lambda s: pd.Series({
            "COUNT": len(s.Sample),
            "PERCENTAGE": f"{len(s.Sample) / n_total_samples:.2%}",
            "SAMPLES": list(s.Sample)
        }), include_groups = False).reset_index().reset_index()

    if CALL_ID in ["CRT", "GCH1"]:
        groupby = Pf8_calls.groupby([f"Pf8: {CALL_ID}"]).apply(lambda s: pd.Series({
            "COUNT": len(s.Sample),
            "PERCENTAGE": f"{len(s.Sample) / n_total_samples:.2%}",
            "SAMPLES": list(s.Sample)
        }), include_groups = False).reset_index().reset_index()

    with sidebar.expander("Show sample call distribution - click on the index to select samples", expanded = True):
        if st.toggle("Set threshold on sample count for each GT combination"):
            gt_combination_sample_count_threshold_input = st.number_input(label = "", value = 100, label_visibility = "collapsed", min_value = 1, max_value = 300)
            groupby = groupby.loc[groupby.COUNT <= gt_combination_sample_count_threshold_input]

            if CALL_ID in ["MDR1", "HRP2", "HRP3", "PM2_PM3"]:
                groupby = groupby.loc[groupby[f"Pf7: {CALL_ID}"] != groupby[f"Pf8: {CALL_ID}"]]

        selected_rows = st.dataframe(
            groupby.drop(columns = ["index", "SAMPLES"]),
            hide_index = True, use_container_width = True, on_select = "rerun", selection_mode = "single-row",
            key = "modify_draft_calls"
        )["selection"]["rows"]
    
    if len(selected_rows) == 0:
        st.stop()

    RAW_INPUT_STR = str(groupby.loc[groupby.index[selected_rows[0]], "SAMPLES"])
    
    sidebar.divider()




if browser_mode == "Inspect final calls":
    # Allow user to choose call region, then display confusion matrix dataframe
    # When the user selects a row, show all samples belonging to that genotype combination
    # e.g., all samples that were GT=1 in Pf7 but GT=0 in Pf8. 
    if "final_coverage_calls.tsv" not in os.listdir("app_files"):
        st.warning("`final_coverage_calls.tsv` not generated. ")
        st.stop()

    n_total_samples = len(Pf8_calls)
    
    calls_comparison = Pf8_calls[["Sample"] + [col for col in Pf8_calls.columns if "Pf7" in col]].rename(columns = {"Sample": "SAMPLE"})
    calls_comparison[["Pf7: CRT", "Pf7: GCH1"]] = "no data"

    calls_comparison = calls_comparison.merge(pd.read_csv("app_files/final_coverage_calls.tsv", sep = "\t"), on = "SAMPLE")

    CALL_IDS = os.listdir(DIAGNOSTIC_PLOT_PATH)
    CALL_ID  = sidebar.selectbox("Call ID:", CALL_IDS)

    groupby = calls_comparison.groupby([f"Pf7: {CALL_ID}", f"{CALL_ID}"]).apply(lambda s: pd.Series({
        "COUNT": len(s.SAMPLE),
        "PERCENTAGE": f"{len(s.SAMPLE) / n_total_samples:.2%}",
        "SAMPLES": list(s.SAMPLE)
    }), include_groups = False).reset_index()

    with sidebar.expander("Show sample call distribution - click on the index to select samples", expanded = True):
        selected_rows = st.dataframe(
            groupby,
            hide_index = True, use_container_width = True, on_select = "rerun", selection_mode = "single-row"
        )["selection"]["rows"]
    
    if len(selected_rows) == 0:
        st.stop()
    
    RAW_INPUT_STR = str(groupby.loc[selected_rows[0], "SAMPLES"])
    
    sidebar.divider()



if browser_mode == "Modification stats":
    # Overkill code to generate an interactive sunburst plot describing how the manual
    # curation process went. 
    local_directory_files = os.listdir("app_files")

    if "draft_coverage_calls.tsv" not in local_directory_files:
        st.warning("`draft_coverage_calls.tsv` not generated. ")
        st.stop()

    if "final_coverage_calls.tsv" not in local_directory_files:
        st.warning("`final_coverage_calls.tsv` not generated. ")
        st.stop()
    
    CALL_IDS = os.listdir(DIAGNOSTIC_PLOT_PATH)
    
    gt_combination_sample_count_threshold_input = sidebar.number_input(label = "GT combination sample count threshold used", value = 100, min_value = 1, max_value = 300)

    pf7_call_column_map = {col: col.partition(": ")[2] for col in Pf8_calls.columns if "Pf7" in col}
    pf7_calls = Pf8_calls[["Sample"] + list(pf7_call_column_map.keys())].rename(columns = pf7_call_column_map).rename(columns = {"Sample": "SAMPLE"})
    pf7_calls_long = pf7_calls.melt(id_vars = ["SAMPLE"], var_name = "CALL_ID", value_name = "PF7_GT")

    draft_calls = pd.read_csv("app_files/draft_coverage_calls.tsv", sep = "\t").set_index("SAMPLE")
    draft_calls_long = draft_calls.reset_index().melt(id_vars = ["SAMPLE"], var_name = "CALL_ID", value_name = "DRAFT_GT")

    final_calls = pd.read_csv("app_files/final_coverage_calls.tsv", sep = "\t").set_index("SAMPLE")
    final_calls_long = final_calls.reset_index().melt(id_vars = ["SAMPLE"], var_name = "CALL_ID", value_name = "FINAL_GT")

    all_calls = (
        pf7_calls_long
        .merge(draft_calls_long, on = ["SAMPLE", "CALL_ID"], how = "outer")
        .merge(final_calls_long, on = ["SAMPLE", "CALL_ID"], how = "outer")
        .fillna("no data")
    )

    all_calls.loc[all_calls.DRAFT_GT != all_calls.FINAL_GT, "MODIFICATION"] = all_calls.loc[all_calls.DRAFT_GT != all_calls.FINAL_GT, "FINAL_GT"]
    all_calls.MODIFICATION = all_calls.MODIFICATION.fillna("NO CHANGE")

    stats_df = all_calls.groupby(["CALL_ID", "PF7_GT", "DRAFT_GT", "MODIFICATION"]).apply(lambda s: pd.Series({
        "COUNT"        : s.SAMPLE.nunique(),
        "MANUALLY_INSPECTED": s.SAMPLE.nunique() <= gt_combination_sample_count_threshold_input
    })).reset_index()

    ids, parents, labels, values, colours, customdata = [], [], [], [], [], []

    genotype_colour_map = {
        0 : "rgba(127,159,222,0.6)",
        1 : "rgba(255,0,0,1)",
        -1: "rgba(0,0,0,0.1)"
    }
    inspection_colour_map = {
        True : "rgba(255,140,39,1)",
        False: "rgba(0,0,0,0.1)"
    }

    # Level 1: CALL_ID
    level1 = stats_df.groupby(["CALL_ID"]).COUNT.sum().reset_index()

    parents += [""] * len(level1)
    ids     += level1.CALL_ID.tolist()
    labels  += level1.CALL_ID.tolist()
    values  += level1.COUNT.tolist()
    colours += ["rgba(0,0,0,0.3)"] * len(level1)
    customdata += level1.COUNT.apply(lambda x: f"Number of samples: {x}").tolist()

    # Level 2: CALL_ID + DRAFT_GT
    level2 = stats_df.groupby(["CALL_ID", "DRAFT_GT"]).COUNT.sum().reset_index()
    level2["CHILD_ID"] = level2.CALL_ID + level2.DRAFT_GT.astype(str)
    level2["COUNT_SUM"] = level2.groupby("CALL_ID").COUNT.transform("sum")
    level2["RATIO"] = level2.COUNT / level2.COUNT_SUM

    level2["TRANSFORMED_COUNT"] = level2.COUNT.apply(lambda x: np.power(x, 0.7))
    level2["TRANSFORMED_COUNT_SUM"] = level2.groupby("CALL_ID").TRANSFORMED_COUNT.transform("sum")
    level2["SCALED_COUNT_LVL2"] = level2.COUNT_SUM * (level2.TRANSFORMED_COUNT / level2.TRANSFORMED_COUNT_SUM)

    parents += level2.CALL_ID.tolist()
    ids     += level2.CHILD_ID.tolist()
    labels  += level2.DRAFT_GT.tolist()
    values  += level2.SCALED_COUNT_LVL2.tolist()
    colours += level2.DRAFT_GT.map(genotype_colour_map).tolist()
    customdata += level2[["CALL_ID", "DRAFT_GT", "COUNT", "COUNT_SUM", "RATIO"]].apply(
        lambda x: f"Proportion of {x[0]} draft calls with GT = {x[1]}: {x[2]} out of {x[3]} ({x[4]:.02%})", axis = 1
    ).tolist()
    
    # Level 3: Manually Inspected
    level3 = stats_df.groupby(["CALL_ID", "DRAFT_GT", "MANUALLY_INSPECTED"]).COUNT.sum().reset_index()
    level3["PARENT_ID"] = level3.CALL_ID + level3.DRAFT_GT.astype(str)
    level3["CHILD_ID"] = level3.PARENT_ID + level3.MANUALLY_INSPECTED.map({True: "Inspected", False: "Not inspected"})
    level3["COUNT_SUM"] = level3.groupby(["CALL_ID", "DRAFT_GT"]).COUNT.transform("sum")
    level3["RATIO"] = level3.COUNT / level3.COUNT_SUM

    level3 = level3.merge(
        level2[["CALL_ID", "DRAFT_GT", "SCALED_COUNT_LVL2"]],
        on = ["CALL_ID", "DRAFT_GT"], how = "outer")
    
    level3["TRANSFORMED_COUNT"] = level3.COUNT.apply(lambda x: np.power(x, 0.5) + 1)
    level3["TRANSFORMED_COUNT_SUM"] = level3.groupby(["CALL_ID", "DRAFT_GT"]).TRANSFORMED_COUNT.transform("sum")
    level3["SCALED_COUNT_LVL3"] = level3.SCALED_COUNT_LVL2 * (level3.TRANSFORMED_COUNT / level3.TRANSFORMED_COUNT_SUM)
    # level3 = level3[level3.MANUALLY_INSPECTED == "Manually inspected"]

    parents += level3.PARENT_ID.tolist()
    ids     += level3.CHILD_ID.tolist()
    labels  += level3.MANUALLY_INSPECTED.map({True: "Inspected", False: "Not inspected"}).tolist()
    values  += level3.SCALED_COUNT_LVL3.tolist()
    colours += level3.MANUALLY_INSPECTED.map(inspection_colour_map).tolist()
    customdata += level3[["CALL_ID", "DRAFT_GT", "MANUALLY_INSPECTED", "COUNT", "COUNT_SUM", "RATIO"]].apply(
        lambda x: f"""Proportion of {x[0]} draft calls with GT= {x[1]} which were {("manually inspected" if x[2] else "not inspected")}: {x[3]} out of {x[4]} ({x[5]:.02%})""", axis = 1
    ).tolist()

    # Level 4: Modifications
    level4 = stats_df.groupby(["CALL_ID", "DRAFT_GT", "MANUALLY_INSPECTED", "MODIFICATION"]).COUNT.sum().reset_index()
    level4["PARENT_ID"] = level4.CALL_ID + level4.DRAFT_GT.astype(str) + level4.MANUALLY_INSPECTED.map({True: "Inspected", False: "Not inspected"})
    level4["CHILD_ID"] = level4.PARENT_ID + level4.MODIFICATION.astype(str)
    level4["COUNT_SUM"] = level4.groupby(["CALL_ID", "DRAFT_GT", "MANUALLY_INSPECTED"]).COUNT.transform("sum")
    level4["RATIO"] = level4.COUNT / level4.COUNT_SUM

    level4 = level4.merge(
        level3[["CALL_ID", "DRAFT_GT", "SCALED_COUNT_LVL3", "MANUALLY_INSPECTED"]],
        on = ["CALL_ID", "DRAFT_GT", "MANUALLY_INSPECTED"], how = "outer")
    
    level4["TRANSFORMED_COUNT"] = level4.COUNT.apply(lambda x: np.power(x, 0.5) + 1)
    level4["TRANSFORMED_COUNT_SUM"] = level4.groupby(["CALL_ID", "DRAFT_GT", "MANUALLY_INSPECTED"]).TRANSFORMED_COUNT.transform("sum")
    level4["SCALED_COUNT_LVL4"] = level4.SCALED_COUNT_LVL3 * (level4.TRANSFORMED_COUNT / level4.TRANSFORMED_COUNT_SUM)
    level4 = level4[level4.MODIFICATION != "NO CHANGE"]

    parents += level4.PARENT_ID.tolist()
    ids     += level4.CHILD_ID.tolist()
    labels  += level4.MODIFICATION.tolist()
    values  += level4.SCALED_COUNT_LVL4.tolist()
    colours += level4.MODIFICATION.map(genotype_colour_map).tolist()
    customdata += level4[["CALL_ID", "DRAFT_GT", "MODIFICATION", "COUNT", "COUNT_SUM", "RATIO"]].apply(
        lambda x: f"""Proportion of manually inspected {x[0]} draft calls with GT= {x[1]}
which were modified to {int(x[2])}: {x[3]} out of {x[4]} ({x[5]:.02%})""", axis = 1
    ).tolist()

    # Create Sunburst plot
    fig = go.Figure(go.Sunburst(
        ids=ids,
        parents=parents,
        labels=labels,
        values=values,
        branchvalues="total",
        hovertemplate="<b>%{label}</b><br>%{customdata}<extra></extra>",
        customdata=customdata
    ))

    fig.update_traces(marker=dict(colors=colours, line=dict(color="black", width=0.2)),
        insidetextorientation="radial", textfont=dict(color="black", family="Arial", size=15)
    )
    fig.update_layout(
        width=1300, height=1100,
        uniformtext_minsize=10, uniformtext_mode=False
    )

    fig.write_html("app_files/final_coverage_results.html")

    main.plotly_chart(fig, use_container_width = True)



if browser_mode == "Inspect modifications":

    local_directory_files = os.listdir("app_files")

    if "draft_coverage_calls.tsv" not in local_directory_files:
        st.warning("`draft_coverage_calls.tsv` not generated. ")
        st.stop()

    if "final_coverage_calls.tsv" not in local_directory_files:
        st.warning("`final_coverage_calls.tsv` not generated. ")
        st.stop()

    CALL_IDS = os.listdir(DIAGNOSTIC_PLOT_PATH)
    CALL_ID  = sidebar.selectbox("Call ID:", CALL_IDS)

    draft_calls = pd.read_csv("app_files/draft_coverage_calls.tsv", sep = "\t").set_index("SAMPLE")
    draft_calls.columns = ["draft:"+col for col in draft_calls.columns]
    final_calls = pd.read_csv("app_files/final_coverage_calls.tsv", sep = "\t").set_index("SAMPLE")
    final_calls.columns = ["final:"+col for col in final_calls.columns]

    merged_calls = draft_calls.join(final_calls).reset_index()

    RAW_INPUT_STR = str(merged_calls.loc[merged_calls[f"draft:{CALL_ID}"] != merged_calls[f"final:{CALL_ID}"]].SAMPLE.tolist())

    sidebar.divider()




if browser_mode in ["General browsing", "Modify draft calls", "Inspect modifications", "Inspect final calls"]:
    # ----- Cache load useful files -----
    plot_paths              = cache_load_plot_paths(CALL_ID, RAW_INPUT_STR)
    failed_ploidy           = cache_load_contig_ploidy_fails()
    goi, chrom_lengths_dict = cache_load_bounds_assets()

    CHROMOSOME = goi.loc[goi.CALL_ID == CALL_ID, "CONTIG"].values[0]

    # ----- Initialise sample indexing to assist with button navigation -----
    SAMPLES              = list(plot_paths.keys())
    SAMPLE_TO_INDEX_DICT = {sample: i for i, sample in enumerate(SAMPLES)}

    if len(SAMPLES) == 0:
        st.warning("No samples found. ")
        st.stop()

    if "sample_index" not in st.session_state:
        st.session_state.sample_index = 0

    # ----- Initialise button navigation and update session state -----
    left, right = sidebar.columns(2)
    if left.button("Previous", use_container_width = True):
        st.session_state.sample_index -= 1

    if right.button("Next", use_container_width = True):
        st.session_state.sample_index += 1

    st.session_state.sample_index = st.session_state.sample_index % len(SAMPLES)

    if st.session_state.sample_index == 0:
        st.toast("⚠️ At sample index 0 ⚠️")

    # ----- Sample selection box -----
    SAMPLE = sidebar.selectbox(
        "Sample:", SAMPLES, index = st.session_state.sample_index,
        on_change = update_sample_index, key = "sample_selectbox")

    # ----- Navigation stats and contig ploidy fail error  -----
    sidebar.html(f"""<h6 style='text-align: center'>Sample index: {st.session_state.sample_index}<br>({len(SAMPLES)} total samples)</h6>""")

    PLOIDY = failed_ploidy.loc[(failed_ploidy.SAMPLE_ID == SAMPLE) & (failed_ploidy.CHROM == CHROMOSOME), "PLOIDY"].values

    if len(PLOIDY) > 0:
        sidebar.html(f"""<h6 style='text-align: center; color: red'>WARNING: Dodgy contig ploidy number: {PLOIDY[0]}</h6>""")

    # ===== Display plot =====
    left_gap, plot, right_panel, right_gap = main.columns([.1, 2, 2, .1])

    SAMPLE = SAMPLES[st.session_state.sample_index]
    plot.image(plot_paths[SAMPLE], use_column_width = True)

    gene_start, gene_end, lower_bound, upper_bound = _get_bounds(CALL_ID, CHROMOSOME)

    if CALL_ID != "CRT":
        sample_broad_cn_ratio_path = os.path.join(BROAD_CN_RATIO_PATH, f"{SAMPLE}.dCR.tsv")

        if os.path.exists(sample_broad_cn_ratio_path):
            broad_dCR = pd.read_csv(sample_broad_cn_ratio_path, sep = "\t", skiprows = 18)
            broad_dCR["INTERVAL_MIDPOINT"] = broad_dCR[["START", "END"]].mean(axis = 1)
            broad_dCR = broad_dCR.groupby("CONTIG").apply(lambda s: {
                "INTERVAL_MIDPOINT": np.array(s.INTERVAL_MIDPOINT),
                "LINEAR_COPY_RATIO": np.array(s.LINEAR_COPY_RATIO)
            }).to_dict()[CHROMOSOME]
        
            fig, ax = plt.subplots(figsize = (7, 2))
            ax.scatter(broad_dCR["INTERVAL_MIDPOINT"],
                    broad_dCR["LINEAR_COPY_RATIO"],
                    c = "r", s = 2, alpha = 0.7, zorder = 1)
            
            ax.set_yticks([1, 2, 3, 4, 5])
            ax.set_xlim(lower_bound, upper_bound)
            ax.set_ylim(-0.05, 5.1)
            ax.set_ylabel("COPY\nNUMBER\nRATIO")
            ax.axvline(x = gene_start, c = "red", alpha = 0.3)
            ax.axvline(x = gene_end,   c = "red", alpha = 0.3)
        
            plt.grid(True, axis = "y")
            plt.tight_layout()
        
            plot.pyplot(fig, use_container_width = True)

    # ===== Calls comparison =====
    calls_dict = Pf7_resistance_calls.loc[
        Pf7_resistance_calls.index == SAMPLE,
        [col for col in Pf7_resistance_calls.columns if CALL_ID in col]
    ].to_dict("records")

    if len(calls_dict) != 0:
        if browser_mode != "Inspect modifications":
            for call_name, call in calls_dict[0].items():
                if call == "":
                    continue

                text_colour = {
                    -1: "gray",
                    0 : "black",
                    1 : "red"
                }[call]
                
                sidebar.html(f"""<h6 style='text-align: center; color: {text_colour}'>{call_name}: {int(call)}</h6>""")
        
        sidebar.divider()


if browser_mode == "Modify draft calls":
    # Allows user to manually modify calls. 
    if "manual changes" not in st.session_state or right_panel.button("Reset manual changes"):
        st.session_state["manual changes"] = {call_id: {} for call_id in CALL_IDS}

    right_panel.divider()
    
    DRAFT_CALL = Pf8_calls.loc[Pf8_calls.Sample == SAMPLE, f"Pf8: {CALL_ID}"].values[0]
    TRACKING_CALL = DRAFT_CALL

    if SAMPLE in st.session_state["manual changes"][CALL_ID]:
        MODIFIED_CALL = st.session_state["manual changes"][CALL_ID][SAMPLE]
        TRACKING_CALL = MODIFIED_CALL
    
    possible_genotypes = [-1, 0, 1]
    FINAL_CALL = right_panel.radio(
        "Change call manually",
        possible_genotypes,
        index = possible_genotypes.index(TRACKING_CALL), # Needed to handle st.session_state properly when revisiting modified calls
        key = f"{CALL_ID}-{SAMPLE}" # Unique key value is needed here to update the position of the index when moving between samples
    )

    if FINAL_CALL != DRAFT_CALL:
        right_panel.html(f"""<h6 style='text-align: center; color: red'>CNV genotype has been changed<br>from {DRAFT_CALL} to {FINAL_CALL}</h6>""")
        st.session_state["manual changes"][CALL_ID][SAMPLE] = FINAL_CALL

    if FINAL_CALL == DRAFT_CALL and SAMPLE in st.session_state["manual changes"][CALL_ID]:
        st.session_state["manual changes"][CALL_ID].pop(SAMPLE)

    right_panel.divider()
    
    if right_panel.toggle("See all changes made so far"):
        right_panel.write(f"""{sum([len(v) for v in st.session_state["manual changes"].values()])} changes so far.""")
        right_panel.write(st.session_state["manual changes"])

    right_panel.divider()

    if right_panel.button("Save `changes_checkpoint.json`"):
        with open("changes_checkpoint.json", "w") as file:
            json.dump(st.session_state["manual changes"], file, indent = 4)
    
    if right_panel.button("Load `changes_checkpoint.json`"):
        with open("changes_checkpoint.json", "r") as file:
            st.session_state["manual changes"] = json.load(file)
    
    if right_panel.button("Convert manual changes to `final_coverage_calls.tsv`"):
        draft_calls_df = Pf8_calls.set_index("Sample")[[f"Pf8: {col}" for col in CALL_IDS]].rename(columns = {f"Pf8: {col}": col for col in CALL_IDS})
        final_calls_df = draft_calls_df.copy(deep = True)

        for call_id in CALL_IDS:
            for sample_to_modify, new_call in st.session_state["manual changes"][call_id].items():
                final_calls_df.loc[sample_to_modify, call_id] = new_call

        n_changes_detected = (final_calls_df != draft_calls_df).sum().sum()
        n_changes_in_dict  = sum([len(st.session_state["manual changes"][call_id]) for call_id in CALL_IDS])

        final_calls_df.reset_index().rename(columns = {"Sample": "SAMPLE"}).to_csv("app_files/final_coverage_calls.tsv", sep = "\t", index = False)
        
        right_panel.success(f"`final_coverage_calls.tsv` created! {n_changes_in_dict} manual changes were made and {n_changes_detected} differences were detected between the dataframes. ")

if browser_mode == "Inspect modifications":
    gt_1 = merged_calls.loc[merged_calls.SAMPLE == SAMPLE, f"draft:{CALL_ID}"].values[0]
    gt_2 = merged_calls.loc[merged_calls.SAMPLE == SAMPLE, f"final:{CALL_ID}"].values[0]
    
    right_panel.html(f"""<h4 style='text-align: center'>Genotype was manually changed<br>from {gt_1} to {gt_2}</h4>""")

    right_panel.divider()