import glob
import os
import random
import re
import time
import pandas as pd
import numpy as np
from Utils import shell_command
import itertools
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from graph_metagenome_v_metagenome import concat_graph
from collections import Counter

def pair_generator(items):
    aps = []
    for i in items:
        for i2 in items:
            p = f"{i}:{i2}"
            if f"{i2}:{i}" not in aps and p not in aps:
                aps.append(p)
    return aps

def qual(quals): # convert ASCII to int Phred score
    phreds = [(ord(char)-33) for char in list(str(quals))]
    return phreds

def indel_count(end_str): # tells indel_prune how many characters to remove for indels
    count = ""
    for j in end_str:
        if j.isnumeric():
            count += j
        else:
            return int(count), len(count)

def indel_prune(read, rquals, mquals, rqual, mqual): # read given as str, removes indels & skips, thresholds q
    refined = list(re.sub("[$/#^]", '', re.sub("(?<=\^).", "", read)))  # delete ^], non-base chars
    rquals, mquals, to_del = list(rquals), list(mquals), []
    for i in range(len(refined)):
        if refined[i] == "-" or refined[i] == "+":
            indel_range, indel_digit = indel_count(refined[(i+1):])
            to_del += list(range(i, i + 1 + indel_digit + indel_range))
    for ele in sorted(list(set(to_del)), reverse=True):
        del refined[ele] # len(read) should match len(qualities) at this point
    if len(set([len(refined), len(rquals), len(mquals)])) > 1:
        print(read, refined, rquals, mquals, len(read), len(refined), len(rquals), len(mquals), flush=True)
        raise KeyError("read, sequencing qualities, and mapping qualities are not same lengths")
    zipped = zip(refined, rquals, mquals)
    zipped = [z for z in zipped if z[0] not in ["<", ">", "*"] and z[1]>rqual and z[2]>mqual]
    if len(zipped)==0:
        return np.nan
    else:
        unzipped = list(zip(*zipped))
        return list(unzipped[0])

def compute(reads_1, reads_2): # takes processed reads from both mpiles, determines if match
    unique_1, unique_2 = set(reads_1), set(reads_2)  # only unique bases in reads from mpile 1
    if unique_1 == unique_2: # at locus, both mpiles have same set of bases
        return "strict_concur"
    elif unique_1.issubset(unique_2) or unique_2.issubset(unique_1):  # one mpile has subset of other
        return "lax_concur"
    else:  # inconsistency between mpiles
        return "mismatch"

def prune(pair, mp1, mp2):
    start_time = time.time()
    print(f"starting comparison for {pair}", flush=True)
    pair = pair.split(":")
    merged = mp1.merge(mp2, on="Merge", how="inner")  # merge pair together so we can compare
    if merged.empty:
        print("no positions shared by mpileups", flush=True)
        return
    elif pair[0]==pair[1]: # compare mpileup to itself
        merged['Label'] = merged.apply(lambda a: compute(a[f"readbases_{pair[0]}_x"], a[f"readbases_{pair[1]}_y"]), axis=1)
    else:
        merged['Label'] = merged.apply(lambda a: compute(a[f"readbases_{pair[0]}"], a[f"readbases_{pair[1]}"]), axis=1)  # returns label of SNP, "concur", "noseq", or "no_VCF"
    with_data = merged.shape[0]
    strict_concur_count = (merged["Label"] == "strict_concur").sum()
    lax_concur_count = (merged["Label"] == "lax_concur").sum() + strict_concur_count
    prune_stats = pd.DataFrame(data={"sample_1": [pair[0]], "sample_2": [pair[1]], "with_data": [with_data],
                                     "strict_concur": [strict_concur_count], "lax_concur": [lax_concur_count],
                                     "strict_score": [(strict_concur_count / with_data) * 100],
                                     "lax_score": [(lax_concur_count / with_data) * 100]})
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = elapsed_time // 60, elapsed_time % 60
    print(f"{pair[0]}:{pair[1]} comparison stats saved; took {minutes} minutes and {seconds:.2f} seconds \n", flush=True)
    return prune_stats

def prep(name, source, tmp, rqual, mqual, ji):
    if os.path.exists(f"{tmp}/tmp_mp/{name}.df.gz") or os.path.exists(f"{tmp}/tmp_mp/{name}.txt"):
        print(f"prep already done {name}", flush=True)
        return
    prep_start = time.time()
    rsync_mp_cmd = f"rsync -avz {source}{name}.pileup.gz {tmp}/{ji}_{name}.pileup.gz"
    shell_command(rsync_mp_cmd)
    print(f"done: {rsync_mp_cmd}")
    mpdf = pd.read_csv(f"{tmp}/{ji}_{name}.pileup.gz", sep="\t")
    mpdf =mpdf[mpdf.apply(lambda row: row["rqual"]==row["rqual"] and row["mqual"]==row["mqual"] and row["rqual"]!="no_seq", axis=1)]
    if mpdf.shape[0] == 0: # no sequencing left
        print("empty, no sequencing after filtering for nans in quality scores", flush=True)
        shell_command(f"echo empty, no sequencing after filtering for nans in quality scores > {tmp}/tmp_mp/{name}.txt")
        shell_command(f"rm {tmp}/{ji}_{name}.pileup.gz")
        return
    mpdf['rqual'] = mpdf['rqual'].apply(lambda q: qual(q))
    mpdf['mqual'] = mpdf['mqual'].apply(lambda w: qual(w))
    print("phreds calculated", flush=True)
    mpdf["readbases"] = mpdf.apply(
        lambda y: y["readbases"].replace(".", y["reference"]).replace(",", y["reference"]),
        axis=1)  # replace . and , with reference
    print(". and , replaced", flush=True)
    mpdf[f"readbases_{name}"] = mpdf.apply(lambda z:
                                        indel_prune(z["readbases"].lower(), z["rqual"], z["mqual"], rqual, mqual),
                                        axis=1)  # decode readbases, remove indels, remove bases w quality < threshold
    mpdf = mpdf[mpdf[f"readbases_{name}"].apply(lambda row: row == row)] # remove loci where all readbases failed
    print("mpileup readbases filtered and converted", flush=True)
    mpdf = mpdf[["Merge",f"readbases_{name}"]] #keep only locus and readbases at locus
    if mpdf.empty:
        print("nothing after second round of filtering", flush=True)
        shell_command(f"echo empty, nothing after second round of filtering > {tmp}/tmp_mp/{name}.txt")
    else:
        mpdf.to_pickle(f"{tmp}/{ji}_{name}_final.df.gz", compression="gzip")
        shell_command(f"rsync -avz {tmp}/{ji}_{name}_final.df.gz {tmp}/tmp_mp/{name}.df.gz")
        print(f"saved: rsync -avz {tmp}/{ji}_{name}_final.df.gz {tmp}/tmp_mp/{name}.df.gz", flush=True)
        shell_command(f"rm {tmp}/{ji}_{name}_final.df.gz")
    shell_command(f"rm {tmp}/{ji}_{name}.pileup.gz")
    prep_end = time.time()
    prep_elapsed = prep_end - prep_start
    prep_minutes, prep_seconds = prep_elapsed // 60, prep_elapsed % 60
    print(f"prep completed for {name}; took {prep_minutes} minutes and {prep_seconds:.2f} seconds \n", flush=True)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: <exe> <path to config file>', flush=True)
        print(len(sys.argv), sys.argv, flush=True)
        sys.exit(-1)
    else:
        _, config_path = sys.argv
        config = pd.read_csv(config_path, sep=":", header=None)
        args = config[1].to_list()
        paths, tmp_path, dest_path, pairs, rqualstr, mqualstr, job_id, graph_option = args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7]
        print("user input from config file: ", paths, tmp_path, dest_path, pairs, rqualstr, mqualstr, job_id, flush=True)
        rquality, mquality, source_path, dest_path = int(rqualstr), int(mqualstr), paths.split("*")[0], f"{dest_path}/SNPsleuth_{job_id}"
        if pairs in ["all","same"] or pairs.endswith(".fam"):
            rei = "^" + source_path + "(.*)" + paths.split("*")[1] + "$"
            samples = [re.match(fr"{rei}", x).group(1) for x in glob.glob(paths)]
            if pairs == ["same"] or pairs.endswith(".fam"):
                indivs = set(["_".join(s.split("_")[:-1]) for s in samples])
                if pairs.endswith(".fam"):
                    famdf = pd.read_csv(pairs, sep="\s+", header=None)
                    famlist = famdf[0].to_list()
                    indivs = [i for i in indivs if i not in famlist]
                    pairs = "same" # functionally the same as "same" option now
                by_indiv = [[s for s in samples if "_".join(s.split("_")[:-1]) == i] for i in indivs]
                samples = []
                for b in by_indiv:
                    if len(b) > 1:
                        samples.extend(b) # only include samples from individuals with 2+ samples
        else: # specific file given with pairs to compare
            pairs_df = pd.read_csv(pairs, sep="\s+", header=None)
            samples = set(pairs_df[0].to_list() + pairs_df[1].to_list())

        # make sure files haven't already been preprocessed
        to_preprocess = [s for s in samples if os.path.exists(f"{tmp_path}/tmp_mp/{s}.df.gz")==False and os.path.exists(f"{tmp_path}/tmp_mp/{s}.txt")==False]
        print(rquality, mquality, source_path, dest_path, samples, len(samples), len(to_preprocess), flush=True)
        # make folders necessary. if folders already exist, will keep them
        to_make = [dest_path, tmp_path, f"{tmp_path}/SNPstats_{job_id}", f"{tmp_path}/tmp_mp"]
        for m in to_make:
            try:
                shell_command(f"mkdir {m}")  # make main folder
                print(f"folder created: {m}", flush=True)
            except:
                print(f"folder already exists: {m}", flush=True)

        # send off preprocessing jobs
        # PARALLELIZE IF DESIRED
        for tp in to_preprocess:
            prep(tp, source_path, tmp_path, rquality, mquality, job_id)
        completed = [g.split("/")[-1].split(".")[0] for g in glob.glob(f"{tmp_path}/tmp_mp/*.df.gz")]
        samples = sorted([smp for smp in samples if smp in completed])
        print(f"number of samples that finished preprocessing: {len(samples)}", flush=True)

        # preprocessing completed, begin comparisons
        if pairs=="all": # all possible pairs of samples given a folder
            all_pairs = pair_generator(samples)
        elif pairs =="same": # compare only samples from the same donors
            combos = [list(itertools.combinations(b, 2)) for b in by_indiv]
            all_pairs = []
            for c in combos:
                all_pairs.extend(c)
            all_pairs = [f"{a[0]}:{a[1]}" for a in all_pairs]
        else: # path to space-delimited file with all desired pairs for comparison, each row has one pair
            pairs_df = pairs_df[pairs_df.apply(lambda row: row[0] in samples and row[1] in samples, axis=1)]
            pairs_df["pair"] = pairs_df[0] + ":" + pairs_df[1]
            all_pairs = set(pairs_df["pair"].to_list())
        print("number of comparisons to run: ", len(all_pairs), flush=True)
        # load all mpileups ahead of time
        print(f"{len(samples)} mpileups to load", flush=True)
        mpile_df_dict = {m: pd.read_pickle(f"{tmp_path}/tmp_mp/{m}.df.gz") for m in samples}
        # initialize growing dataframe
        concat = pd.DataFrame(columns=["sample_1", "sample_2", "with_data", "strict_concur", "lax_concur", "strict_score", "lax_score"])
        # PARALLELIZE IF DESIRED
        for ap in all_pairs:
            concat = pd.concat([concat, prune(ap, mpile_df_dict[ap.split(":")[0]], mpile_df_dict[ap.split(":")[1]])])

        covg_search = glob.glob(f"{'/'.join(dest_path.split('/')[:-1])}/preprocess*/coverage.tsv")
        if len(covg_search) == 1:
            print(f"coverage file exists: {covg_search[0]}", flush=True)
            shell_command(f"rsync -avz {covg_search[0]} {tmp_path}/SNPstats_{job_id}/")
            covg = pd.read_csv(f"{tmp_path}/SNPstats_{job_id}/coverage.tsv", sep="\t")
            covg.rename(columns={"coverage":"coverage_1","sample":"sample_1"}, inplace=True)
            concat = concat.merge(covg, on="sample_1", how="left")
            covg.rename(columns={"coverage_1":"coverage_2","sample_1":"sample_2"}, inplace=True)
            concat = concat.merge(covg, on="sample_2", how="left")
            covg = pd.read_csv(f"{tmp_path}/SNPstats_{job_id}/coverage.tsv", sep="\t")
        else:
            covg = pd.DataFrame()
        concat.to_csv(f"{tmp_path}/SNPstats_{job_id}/concat_{job_id}.tsv", sep="\t", index=False)
        rsync_cat_cmd = f"rsync -avz {tmp_path}/SNPstats_{job_id}/concat_{job_id}.tsv {dest_path}/"
        print(rsync_cat_cmd, flush=True)
        shell_command(rsync_cat_cmd)

        if graph_option=="yes":
            concat_graph(concat, tmp_path, job_id, covg)
            rsync_final_cmd = f"rsync -avz {tmp_path}/SNPstats_{job_id}/ {dest_path}/"
            print(rsync_final_cmd, flush=True)
            shell_command(rsync_final_cmd)
            print(f"transferred files to final destination", flush=True)

        shell_command(f"rm -rf {tmp_path}/SNPstats_{job_id}")
        shell_command(f"rm -f {tmp_path}/*.*")
        print("removed tmp files", flush=True)
