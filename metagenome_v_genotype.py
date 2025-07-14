import glob
import re
import os
import pandas as pd
import numpy as np
import time
from Utils import shell_command
import itertools
import sys
import time
import matplotlib.pyplot as plt
import seaborn as sns
import random
from graph_metagenome_v_genotype import concat_graph
import gc

def pair_generator(ss, gs):
    aps = []
    # all same-genotype comparisons together
    for g in gs:
        for s in ss:
            p = f"{s}:{g}"
            aps.append(p)
    return aps

def chunker(items, csize):
    # random.shuffle(items) , not randomizing so that genotypes stay together
    tmp_chunk = list(chunked(items, int(csize)))
    return tuple(tuple(l) for l in tmp_chunk)

def qual(quals): # convert ASCII to int Phred score
    phreds = [(ord(char)-33) for char in list(str(quals))]
    return phreds

def genotype(ref, alts, gtype): # converts genotype from format #/# to set of base letters
    if gtype.split("/")[0].isnumeric() == False or gtype.split("/")[1].isnumeric() == False:
        return np.nan
    else:
        choices = [ref] + alts.split(",")  # list of choices, list[0] is reference
        types = [choices[int(base)].lower() for base in gtype.split("/")]
        return set(types)  # list of bases (lowercase) matching genotype

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

def compute(reads, gtype): # takes processed read, determines if match or SNP
    unique = set(reads) # only unique bases in reads
    if unique==gtype:
        return "strict_concur"
    elif unique.issubset(gtype):  # read bases = or are subset of genotype
        return "lax_concur"
    else:  # inconsistent: base in readbases not in genotype
        return "mismatch"

def prune(name, mpile_df, vcf_df):
    start_time = time.time()
    print(f"starting {name} at {start_time}", flush=True)
    pair_mpile, pair_vcf = name.split(":")[0], name.split(":")[1]
    indiv_vcf = pair_vcf.split("_")[0] + "_" + pair_vcf.split("_")[-1]
    merged = mpile_df.merge(vcf_df, on="Merge", how="inner")
    merged["Label"] = merged.apply(lambda row: compute(row[f"readbases_{pair_mpile}"], row[indiv_vcf]), axis=1)
    # save stats df
    with_data = merged.shape[0] # positions with both mpile and VCF
    strict_concur_count = (merged["Label"] == "strict_concur").sum() # positions where mpile and VCF concur
    lax_concur_count = (merged["Label"] == "lax_concur").sum() + strict_concur_count  # positions where mpile and VCF concur
    prune_stats = pd.DataFrame(data={"sample": [pair_mpile], "genotype": [pair_vcf],
                                     "with_data": [with_data], "strict_concur": [strict_concur_count],
                                     "lax_concur": [lax_concur_count],
                                     "strict_score": [(strict_concur_count/with_data)*100],
                                     "lax_score": [(lax_concur_count/with_data)*100]})
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = elapsed_time // 60, elapsed_time % 60
    print(f"{name} stats saved; took {minutes} minutes and {seconds:.2f} seconds \n", flush=True)
    return prune_stats

def vcfprep(name, tmp, vcfp, ji):
    if os.path.exists(f"{tmp}/tmp_vcf/{name}.df.gz"):
        print(f"already done {name}", flush=True)
        return
    vcfprep_start = time.time()
    indiv_plink = name.split("_")[0] + " " + name.split("_")[-1]
    indiv_vcf = name.split("_")[0] + "_" + name.split("_")[-1]
    shell_command(f"echo {indiv_plink} > {tmp}/{ji}_{name}_plinkkeep.txt")
    shell_command(f"plink --bfile {vcfp} --keep {tmp}/{ji}_{name}_plinkkeep.txt --recode vcf --out {tmp}/{ji}_{name}_plink")
    shell_command(f"bcftools view -v snps -Oz -o {tmp}/{ji}_{name}_plink_snps.vcf.gz {tmp}/{ji}_{name}_plink.vcf")
    shell_command(f"bcftools head {tmp}/{ji}_{name}_plink_snps.vcf.gz > {tmp}/{ji}_{name}_plink_snps_head.txt")
    with open(f"{tmp}/{ji}_{name}_plink_snps_head.txt", 'r') as fp:
        lines = len(fp.readlines()) - 1
    df = pd.read_csv(f"{tmp}/{ji}_{name}_plink_snps.vcf.gz", sep="\t", skiprows=lines, usecols=["#CHROM", "POS", "REF", "ALT", indiv_vcf])
    df["Merge"] = df["#CHROM"].astype(str) + "_" + df["POS"].astype(str)
    df[indiv_vcf] = df.apply(lambda row: genotype(row["REF"], row["ALT"], row[indiv_vcf]), axis=1)
    df = df[df[indiv_vcf]==df[indiv_vcf]]
    df = df[["Merge", indiv_vcf]]
    df.to_pickle(f"{tmp}/{ji}_{name}_plink_final.df.gz", compression="gzip")
    shell_command(f"rsync -avz {tmp}/{ji}_{name}_plink_final.df.gz {tmp}/tmp_vcf/{name}.df.gz")
    shell_command(f"rm {tmp}/{ji}_{name}*plink*")
    vcfprep_end = time.time()
    vcfprep_elapsed = vcfprep_end - vcfprep_start
    vcfprep_minutes, vcfprep_seconds = vcfprep_elapsed // 60, vcfprep_elapsed % 60
    print(f"prep completed for {name} vcf; took {vcfprep_minutes} minutes and {vcfprep_seconds:.2f} seconds \n",
          flush=True)

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
        print(config_path, flush=True)
        config = pd.read_csv(config_path, sep=":", header=None)
        args = config[1].to_list()
        paths, tmp_path, dest_path, pairs, vcfpath, rqualstr, mqualstr, job_id, graph_option = args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8]
        print("user input from config file: ", paths, tmp_path, dest_path, pairs, vcfpath, rqualstr, mqualstr, job_id, flush=True)
        pathsnoslash = [i[:-1] if i.endswith("/") else i for i in [tmp_path, dest_path]] # remove end slash if present
        tmp_path, dest_path = pathsnoslash[0], pathsnoslash[1]
        rquality, mquality, source_path, dest_path = int(rqualstr), int(mqualstr), paths.split("*")[0], f"{dest_path}/SNPsleuth_{job_id}"
        fam = pd.read_csv(f"{vcfpath}.fam", sep=" ", header=None)
        fam["fam"] = fam[0]+"_"+fam[1]
        fam = fam["fam"].to_list()
        if pairs in ["self","all","vcf"]:
            rei = "^" + source_path + "(.*)" + paths.split("*")[1] + "$"
            samples = [re.match(fr"{rei}", x).group(1) for x in glob.glob(paths)]
            if pairs=="vcf":
                genos = [f.split("_")[0] for f in fam]
                nogeno = []
            else:
                genos = list(set(["_".join(j.split("_")[:-1]) if len(j.split("_"))>=2 else j for j in samples]))
                if len(genos[0].split("_"))>=2: # include only genotypes that are in the .fam file
                    nogeno = [g for g in genos if g not in fam]
                    genos = [g for g in genos if g in fam]
                else:
                    nogeno = [g for g in genos if f"{g}_{g}" not in fam]
                    genos = [g for g in genos if f"{g}_{g}" in fam]
        else:
            pairs_df = pd.read_csv(pairs, sep="\s+", header=None)
            geno_example = pairs_df.iat[0,1]
            if len(geno_example.split("_"))>=2: # include only pairs with genotypes in the .fam file
                nogeno = pairs_df[pairs_df[1].apply(lambda row: row not in fam)]
                pairs_df = pairs_df[pairs_df[1].apply(lambda row: row in fam)]
            else:
                nogeno = pairs_df[pairs_df[1].apply(lambda row: f"{row}_{row}" not in fam)]
                pairs_df = pairs_df[pairs_df[1].apply(lambda row: f"{row}_{row}" in fam)]
            nogeno = nogeno[1].to_list()
            samples, genos = set(pairs_df[0].to_list()), set(pairs_df[1].to_list())
        print(f"no array genotypes for the following individuals: ", nogeno, flush=True)
        to_preprocess = [s for s in samples if os.path.exists(f"{tmp_path}/tmp_mp/{s}.df.gz") == False and os.path.exists(f"{tmp_path}/tmp_mp/{s}.txt") == False]
        to_preprocess_vcf = [g for g in genos if os.path.exists(f"{tmp_path}/tmp_vcf/{g}.df.gz") == False]
        # make folders necessary. if folders already exist, will keep them
        to_make = [dest_path, tmp_path, f"{tmp_path}/SNPstats_{job_id}", f"{tmp_path}/tmp_mp", f"{tmp_path}/tmp_vcf"]
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

        # send off vcf preprocessing jobs
        # PARALLELIZE IF DESIRED
        for tpv in to_preprocess_vcf:
            vcfprep(tpv, tmp_path, vcfpath, job_id)

        # determine pairs to compare
        if pairs == "self":
            all_pairs = pd.DataFrame(data={"sample":samples})
            all_pairs["genotype"] = all_pairs["sample"].apply(lambda j: '_'.join(j.split('_')[:-1] if len(j.split('_')) >= 2 else f"{j}:{j}")) # find genotype
            all_pairs = all_pairs[all_pairs["genotype"].apply(lambda g: g in genos)] # make sure genotype in fam file
            all_pairs["pair"] = all_pairs["sample"] + ":" + all_pairs["genotype"] # create pair list
            all_pairs.sort_values(by=["genotype"], inplace=True) # sort by genotype
            all_pairs = all_pairs["pair"].to_list() # convert to list
        elif pairs=="all" or pairs=="vcf":  # all
            all_pairs = pair_generator(samples, genos)
        else: # text file with pairs
            pairs_df = pairs_df[pairs_df[0].apply(lambda row: row in samples)] # remove pairs where sample failed filter
            pairs_df.sort_values(by=1, inplace=True) # sort by genotype
            pairs_df["pair"] = pairs_df[0] + ":" + pairs_df[1]
            all_pairs = pairs_df["pair"].to_list()
        print(f"pairs to process: {len(all_pairs)}", flush=True)

        # send off comparison jobs
        # load all mpileups ahead of time
        print(f"{len(samples)} mpileups to load", flush=True)
        mpile_df_dict = {m: pd.read_pickle(f"{tmp_path}/tmp_mp/{m}.df.gz") for m in samples}
        geno_to_load = ""
        # initialize empty dataframe
        concat = pd.DataFrame(columns=["sample_1", "genotype", "with_data", "strict_concur", "lax_concur", "strict_score", "lax_score"])
        # PARALLELIZE IF DESIRED
        for ap in all_pairs:
            if geno_to_load != ap.split(":")[1]: # load in new genotype
                geno_to_load = ap.split(":")[1]
                start_load = time.time()
                geno_loaded = pd.read_pickle(f"{tmp_path}/tmp_vcf/{geno_to_load}.df.gz")
                end_load = time.time()
                elapsed_load = end_load - start_load
                minutes, seconds = elapsed_load // 60, elapsed_load % 60
                print(f"{geno_to_load} loaded; took {minutes} minutes and {seconds:.2f} seconds \n", flush=True)
                gc.collect()
            concat = pd.concat([concat, prune(ap, mpile_df_dict[ap.split(":")[0]], geno_loaded)])
        concat.to_csv(f"{tmp_path}/SNPstats_{job_id}/concat_{job_id}.tsv", sep="\t", index=False)
        covg_search = glob.glob(f"{'/'.join(dest_path.split('/')[:-1])}/preprocess*/coverage.tsv")
        if len(covg_search)==1:
            print(f"coverage file exists: {covg_search[0]}", flush=True)
            shell_command(f"rsync -avz {covg_search[0]} {tmp_path}/SNPstats_{job_id}/")
            covg = pd.read_csv(f"{tmp_path}/SNPstats_{job_id}/coverage.tsv", sep="\t")
            concat = concat.merge(covg, on="sample", how="left")
        else:
            covg = pd.DataFrame()
        concat.to_csv(f"{tmp_path}/SNPstats_{job_id}/concat_{job_id}.tsv", sep="\t", index=False)
        rsync_cat_cmd = f"rsync -avz {tmp_path}/SNPstats_{job_id}/concat_{job_id}.tsv {dest_path}/"
        print(rsync_cat_cmd, flush=True)
        shell_command(rsync_cat_cmd)

        # graph it
        if graph_option == "yes":
            concat_graph(concat, tmp_path, job_id, covg)
            shell_command(f"rsync -avz {tmp_path}/SNPstats_{job_id}/ {dest_path}/")
            print(f"transferred graphs  to final destination", flush=True)

        shell_command(f"rm -rf {tmp_path}/SNPstats_{job_id}")
        shell_command(f"rm -rf {tmp_path}/*.*")
        print(f"removed temporary folder {tmp_path}/SNPstats_{job_id}", flush=True)


