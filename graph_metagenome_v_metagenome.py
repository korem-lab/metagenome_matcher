import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
plt.rcParams.update({'font.size': 15})
import os

def basic_heatmap(matrix, cdf, vm, savepath):
    fig, ((ax_covg_cbar, ax_covg, ax_heatmap, ax_cbar)) = plt.subplots(1, 4, figsize=(24, 20),gridspec_kw={'width_ratios': [1, 1, 11, 1]})
    g = sns.heatmap(matrix, ax=ax_heatmap, cbar_ax=ax_cbar, xticklabels=False, yticklabels=False)
    g.set_facecolor('xkcd:teal')
    sns.heatmap(cdf, xticklabels=False, yticklabels=False, ax=ax_covg, cbar_ax=ax_covg_cbar,
                cmap="viridis", vmax=vm)
    ax_covg_cbar.yaxis.tick_left()
    plt.subplots_adjust(wspace=0.05)
    plt.savefig(savepath)
    plt.clf()

def heatmap_maker(matrix, granges, vlabels, cdf, vm, savepath):
    whole = granges[-1][-1]
    fig, ((ax_covg_cbar, ax_covg, ax_heatmap, ax_cbar),(ax_blank_3, ax_blank_2, ax_legend, ax_blank))= plt.subplots(2,4, figsize=(24, 20), gridspec_kw={'height_ratios': [10, 1],'width_ratios': [1,1,11, 1]})
    g = sns.heatmap(matrix, ax=ax_heatmap, cbar_ax=ax_cbar, xticklabels=False, yticklabels=False)
    g.set_facecolor('xkcd:teal')
    sns.heatmap(cdf, xticklabels=False, yticklabels=False, ax=ax_covg, cbar_ax=ax_covg_cbar,
                cmap="viridis", vmax=vm)
    for i, (start, end) in enumerate(granges):
        ax_legend.axhspan(0, 1, xmin=start / whole, xmax=end / whole, color=f'C{i}', zorder=0)
    xtx = [((x[0] + x[1]) / 2) / whole for x in granges]
    ax_legend.set_yticks([])
    ax_legend.set_yticklabels([])
    ax_legend.set_xticks(xtx)
    ax_legend.set_xticklabels(vlabels)
    ax_covg_cbar.yaxis.tick_left()
    for pos in ['right', 'top', 'bottom', 'left']:
        ax_legend.spines[pos].set_visible(False)
    for x in [ax_blank, ax_blank_2, ax_blank_3]:
        x.set_yticks([])
        x.set_yticklabels([])
        x.set_xticks([])
        x.set_xticklabels([])
        for pos in ['right', 'top', 'bottom', 'left']:
            x.spines[pos].set_visible(False)
    plt.subplots_adjust(hspace=0.05, wspace=0.05)
    plt.savefig(savepath)
    plt.clf()

def matrix_maker(sset, cat): # intakes list of samples and a dataframe with swapped results
    samps_len = len(sset) # how many samples
    strict_matrix, lax_matrix = [0] * samps_len, [0] * samps_len # create empty lists (columns of matrices)
    for s in range(samps_len): # for each sample
        samp = sset[s]
        strict_list, lax_list = [0] * samps_len, [0] * samps_len # create empty lists (rows of matrices)
        slic = cat[cat.apply(lambda row: samp == row["sample_1"] or samp == row["sample_2"], axis=1)] # all rows of df containing sample
        for s2 in range(samps_len): #for each sample
            samp2 = sset[s2]
            if samp2==samp:
                strict_list[s2] = 100
                lax_list[s2] = 100
            else:
                slic2 = slic[slic.apply(lambda row: samp2 == row["sample_2"] or samp2 == row["sample_1"], axis=1)] # all rows of df containing sample 2
                if slic2.empty: # no comparison completed between sample 1 and sample 2
                    strict_list[s2] = np.nan
                    lax_list[s2] = np.nan
                else:
                    strict_list[s2] = slic2["strict_score"].to_list()[0]
                    lax_list[s2] = slic2["lax_score"].to_list()[0]
        strict_matrix[s], lax_matrix[s] = strict_list, lax_list
    strict_matrix = pd.DataFrame(data=np.array(strict_matrix), index=sset, columns=sset)
    lax_matrix = pd.DataFrame(data=np.array(lax_matrix), index=sset, columns=sset)
    return strict_matrix, lax_matrix

def classify_graph(s1, s2): # are two samples from the same individual?
    if s1==s2:
        return "same_sample"
    elif "_".join(s1.split("_")[:-1]) == "_".join(s2.split("_")[:-1]):
        return "same_individual"
    else:
        return "other"

def concat_graph(concat, temp, ji, covg_df):
    samps = list(set(concat["sample_1"].to_list() + concat["sample_2"].to_list()))
    nosame = concat[concat["sample_1"]!=concat["sample_2"]]
    if len(samps[0].split("_"))>=2: # potentially multiple samples from the same individual
        nosame["same"] = nosame.apply(lambda row: classify_graph(row["sample_1"], row["sample_2"]), axis=1)
        # plot boxplots for strict score, differentiating by whether samples came from same vs other person
        plt.figure(figsize=(30, 12))
        sns.boxplot(data=nosame, x="same", y="strict_score", orient= "v", palette='Greys')
        sns.swarmplot(data=nosame, x="same", y="strict_score", orient="v", hue="with_data", palette='rainbow',dodge=False, alpha=0.8)
        plt.savefig(f"{temp}/SNPstats_{ji}/strict_boxplot_samevsother")
        plt.clf()
        # plot boxplots for lax score, differentiating by whether samples came from same vs other person
        plt.figure(figsize=(30, 12))
        sns.boxplot(data=nosame, x="same", y="lax_score", orient= "v", palette='Greys')
        sns.swarmplot(data=nosame, x="same", y="lax_score", orient="v", hue="with_data", palette='rainbow',dodge=False, alpha=0.8)
        plt.savefig(f"{temp}/SNPstats_{ji}/lax_boxplot_samevsother")
        plt.clf()
    # plot boxplot for strict score
    sns.boxplot(data=nosame, x="strict_score", palette="Greys")
    sns.swarmplot(data=nosame, x="strict_score", hue="with_data", palette='rainbow',dodge=False, alpha=0.8)
    plt.savefig(f"{temp}/SNPstats_{ji}/strict_boxplot")
    plt.clf()
    # plot boxplot for lax score
    sns.boxplot(data=nosame, x="lax_score", palette="Greys")
    sns.swarmplot(data=nosame, x="lax_score", hue="with_data", palette='rainbow',dodge=False, alpha=0.8)
    plt.savefig(f"{temp}/SNPstats_{ji}/lax_boxplot")
    plt.clf()
    # make matrix
    concat = concat[concat["with_data"] > 1000]
    samps = list(set(concat["sample_1"].to_list() + concat["sample_2"].to_list()))
    samps.sort()  # sort samples alphabetically
    if covg_df.empty:  # no coverage file provided, instead use size of mpileup as proxy for coverage
        sizes = []
        for s in samps:
            tmp_size_df = pd.read_pickle(f"{temp}/tmp_mp/{s}.df.gz")
            sizes.append(tmp_size_df.shape[0])
        covg_df = pd.DataFrame(data={"sample": samps, "coverage": sizes})
        vmax = 500000
    else:
        vmax = 3
    covg_df.sort_values(by=["sample"], inplace=True)
    covg_df_onecol = covg_df[["coverage"]]
    strict_x, lax_x = matrix_maker(samps, concat) # create matrices for heatmap
    strict_x.to_pickle(f"{temp}/SNPstats_{ji}/strict_matrix.df.gz", compression="gzip")
    lax_x.to_pickle(f"{temp}/SNPstats_{ji}/lax_matrix.df.gz", compression="gzip")
    # graph strict heatmap
    basic_heatmap(strict_x, covg_df_onecol, vmax, f"{temp}/SNPstats_{ji}/strict_heatmap")
    # graph lax heatmap
    basic_heatmap(lax_x, covg_df_onecol, vmax, f"{temp}/SNPstats_{ji}/lax_heatmap")
    # graph heatmap with visit legend
    if len(samps[0].split("_")) >= 2: # we know which individual each sample came from, potentially multiple samples from 1 individual
        indivs = ["_".join(s.split("_")[:-1]) for s in samps]
        indiv_dict = Counter(indivs)  # counter shows how many samples each individual has
        vals = sorted(indiv_dict.values())  # sorted list of number of visits per individual
        tmp_storage, vislabels, group_ranges, range_tracker = [], [], [], 0
        for v in range(vals[-1], 0, -1):  # count backward from greatest number of visits
            tmp = [s for s in samps if indiv_dict["_".join(s.split("_")[:-1])] == v]  # list of samples coming from individual with v samples provided
            if len(tmp) > 0:  # 1+ samples coming from individual with v samples provided
                tmp_storage = tmp_storage + sorted(tmp)  # sort samples alphabetically and add to growing list
                group_ranges.append((range_tracker, (len(tmp) + range_tracker)))
                range_tracker = len(tmp) + range_tracker
                if v != 1:
                    vislabels.append(f"{str(v)} samples")  # label for graph with number of visits
                else:
                    vislabels.append("1 sample")  # label for graph with number of visits
        strict_x, lax_x = matrix_maker(tmp_storage, concat)
        tmp_storage_df = pd.DataFrame(data={"sample": tmp_storage})
        tmp_storage_covg = covg_df.merge(tmp_storage_df, how="right", on="sample")
        tmp_storage_covg = tmp_storage_covg[["coverage"]]
        heatmap_maker(strict_x, group_ranges, vislabels, tmp_storage_covg, vmax, f"{temp}/SNPstats_{ji}/strict_heatmap_visitslabeled")
        heatmap_maker(lax_x, group_ranges, vislabels, tmp_storage_covg, vmax, f"{temp}/SNPstats_{ji}/lax_heatmap_visitslabeled")

