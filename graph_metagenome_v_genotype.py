import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
plt.rcParams.update({'font.size': 15})

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
    fig, ((ax_covg_cbar, ax_covg, ax_heatmap, ax_cbar), (ax_blank_3, ax_blank_2, ax_legend, ax_blank)) = plt.subplots(2,4,figsize=(12,10), gridspec_kw={'height_ratios': [10,1], 'width_ratios': [1,1,11,1]})
    g = sns.heatmap(matrix, xticklabels=False, yticklabels=False, ax=ax_heatmap, cbar_ax=ax_cbar)
    g.set_facecolor('xkcd:teal')
    sns.heatmap(cdf, xticklabels=False, yticklabels=False, ax=ax_covg, cbar_ax=ax_covg_cbar, cmap="viridis",
                vmax=vm)
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

def matrix_maker(sset, gset, cat):
    samps_len, genos_len = len(sset), len(gset)
    strict_matrix, lax_matrix = [0] * samps_len, [0] * samps_len
    for s in range(samps_len):
        samp = sset[s]
        strict_list, lax_list = [0] * genos_len, [0] * genos_len
        slic = cat[cat["sample"].apply(lambda row: samp == row)]
        for g in range(genos_len):
            gt = gset[g]
            slic2 = slic[slic["genotype"].apply(lambda row: gt == row)]
            if slic2.empty:
                strict_list[g] = np.nan
                lax_list[g] = np.nan
            else:
                strict_list[g] = slic2["strict_score"].to_list()[0]
                lax_list[g] = slic2["lax_score"].to_list()[0]
        strict_matrix[s], lax_matrix[s] = strict_list, lax_list
    strict_matrix = pd.DataFrame(data=np.array(strict_matrix), index=sset, columns=gset)
    lax_matrix = pd.DataFrame(data=np.array(lax_matrix), index=sset, columns=gset)
    return strict_matrix, lax_matrix

def concat_graph(concat, temp, ji, covg_df):
    ss, gs = list(concat["sample"].unique()), list(concat["genotype"].unique())
    if len(ss[0].split("_")) >= 2:  # potentially multiple samples from the same individual
        concat["donor"] = concat.apply(lambda row: "donor_individual" if "_".join(row["sample"].split("_")[:-1]) == row["genotype"] else "other_individual", axis=1)
        sns.boxplot(data=concat, x="donor", y="strict_score", orient="v", palette='Greys')
        sns.swarmplot(data=concat, x="donor", y="strict_score", orient="v", hue="with_data", palette='rainbow',dodge=False, alpha=0.8)
        plt.savefig(f"{temp}/SNPstats_{ji}/strict_boxswarm_bydonor")
        plt.clf()
        sns.boxplot(data=concat, x="donor", y="lax_score", orient="v", palette='Greys')
        sns.swarmplot(data=concat, x="donor", y="lax_score", orient="v", hue="with_data", palette='rainbow',dodge=False, alpha=0.8)
        plt.savefig(f"{temp}/SNPstats_{ji}/lax_boxswarm_bydonor")
        plt.clf()
    sns.boxplot(data=concat, x="strict_score", palette='Greys')
    sns.swarmplot(data=concat, x="strict_score", hue="with_data", palette='rainbow',dodge=False, alpha=0.8)
    plt.savefig(f"{temp}/SNPstats_{ji}/strict_boxswarm")
    plt.clf()
    sns.boxplot(data=concat, x="lax_score", palette='Greys')
    sns.swarmplot(data=concat, x="lax_score", hue="with_data", palette='rainbow',dodge=False, alpha=0.8)
    plt.savefig(f"{temp}/SNPstats_{ji}/lax_boxswarm")
    plt.clf()
    # make matrix
    ss.sort()
    gs.sort()
    if covg_df.empty: # no coverage file provided, instead use size of mpileup as proxy for coverage
        sizes = []
        for s in ss:
            tmp_size_df = pd.read_pickle(f"{temp}/tmp_mp/{s}.df.gz")
            sizes.append(tmp_size_df.shape[0])
        covg_df = pd.DataFrame(data={"sample": samps, "coverage": sizes})
        vmax = 500000
    else:
        vmax=3
    covg_df.sort_values(by=["sample"], inplace=True)
    covg_df_onecol = covg_df[["coverage"]]
    strict_x, lax_x = matrix_maker(ss, gs, concat)
    strict_x.to_pickle(f"{temp}/SNPstats_{ji}/strict_matrix.df.gz", compression="gzip")
    lax_x.to_pickle(f"{temp}/SNPstats_{ji}/lax_matrix.df.gz", compression="gzip")
    basic_heatmap(strict_x, covg_df_onecol, vmax, f"{temp}/SNPstats_{ji}/strict_heatmap")
    basic_heatmap(lax_x, covg_df_onecol, vmax, f"{temp}/SNPstats_{ji}/lax_heatmap")
    # only graph samples that have their genotype represented
    ss = [s for s in ss if "_".join(s.split("_")[:-1]) in gs]
    indivs = ["_".join(s.split("_")[:-1]) for s in ss]
    indiv_dict = Counter(indivs)
    vals = sorted(indiv_dict.values())
    tmp_storage, tmp_geno, vislabels, group_ranges, range_tracker = [], [], [], [], 0
    for v in range(vals[-1], 0, -1):
        tmp = [s for s in ss if indiv_dict["_".join(s.split("_")[:-1])] == v]
        if len(tmp) > 0:
            tmp_storage = tmp_storage + sorted(tmp)
            included_genos = list(set(["_".join(t.split("_")[:-1]) for t in tmp]))
            tmp_geno = tmp_geno + sorted(included_genos)
            group_ranges.append((range_tracker, (len(included_genos) + range_tracker)))
            range_tracker = len(included_genos) + range_tracker
            if v != 1:
                vislabels.append(f"{str(v)} samples")
            else:
                vislabels.append("1 sample")
    strict_x, lax_x = matrix_maker(tmp_storage, tmp_geno, concat)
    tmp_storage_df = pd.DataFrame(data={"sample": tmp_storage})
    tmp_storage_covg = covg_df.merge(tmp_storage_df, how="right", on="sample")
    tmp_storage_covg = tmp_storage_covg[["coverage"]]
    strict_x.to_pickle(f"{temp}/SNPstats_{ji}/strict_matrix_withgeno.df.gz", compression="gzip")
    lax_x.to_pickle(f"{temp}/SNPstats_{ji}/lax_matrix_withgeno.df.gz", compression="gzip")
    heatmap_maker(strict_x, group_ranges, vislabels, tmp_storage_covg, vmax, f"{temp}/SNPstats_{ji}/strict_heatmap_withgeno")
    heatmap_maker(lax_x, group_ranges, vislabels, tmp_storage_covg, vmax, f"{temp}/SNPstats_{ji}/lax_heatmap_withgeno")
    
