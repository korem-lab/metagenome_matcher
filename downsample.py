from Utils import shell_command
import os
import re
import glob
import sys
import pandas as pd
import time
import random

def downsample(ogname, fixed, avail_coverage, source, temp, destination_jobID, thresh, seed):
    if os.path.exists(f"{destination_jobID}/seed{seed}/downsample{thresh}/{fixed}.bam"):
        print(f"{destination_jobID}/seed{seed}/downsample{thresh}/{fixed}.bam already done")
        return
    print(f"working on: {destination_jobID}/seed{seed}/downsample{thresh}/{fixed}.bam", flush=True)
    start_time = time.time()
    # move to local/temp space and fix name
    shell_command(f"rsync -avz {source}/{ogname}.bam.gz {temp}/{fixed}_seed{seed}_ds{thresh}_raw.bam.gz")

    # gunzip
    shell_command(f"gunzip {temp}/{fixed}_seed{seed}_ds{thresh}_raw.bam.gz")

    # calculate what fraction of file you should need (we want to chop file as much as possible to speed up subsampling)
    desired_coverage = (int(thresh)*0.06)/1000000 # convert highest desired read threshold to coverage
    fraction = (desired_coverage*1.5)/avail_coverage # what fraction of reads needed to obtain three times the top threshold?
    if fraction < 1: # we have more than 3x the top threshold of reads
        shell_command(f"samtools view -f 3 --subsample {fraction} --subsample-seed {int(seed)} -h -o {temp}/{fixed}_seed{seed}_ds{thresh}_chop.sam {temp}/{fixed}_seed{seed}_ds{thresh}_raw.bam")
    else: # we have 3x the top threshold of reads or less --> do NOT chop, but still filter for perfect pairs
        shell_command(f"samtools view -f 3 -h -o {temp}/{fixed}_seed{seed}_ds{thresh}_chop.sam {temp}/{fixed}_seed{seed}_ds{thresh}_raw.bam")

    # isolate header so that we know how many rows to skip when loading in .sam file as pandas dataframe
    shell_command(f"samtools view -H -o {temp}/{fixed}_seed{seed}_ds{thresh}_chop_header.txt {temp}/{fixed}_seed{seed}_ds{thresh}_chop.sam")
    with open(f"{temp}/{fixed}_seed{seed}_ds{thresh}_chop_header.txt", 'r') as header_file:
        lines = len(header_file.readlines())

    # load in .sam file as dataframe
    sam_df = pd.read_csv(f"{temp}/{fixed}_seed{seed}_ds{thresh}_chop.sam", sep="fakeseparator", skiprows=(lines - 1), header=None, names=["temp"])
    if sam_df.empty:
        print(f"{temp}/{fixed}_seed{seed}_ds{thresh}_chop.sam is truncated", flush=True)
    elif sam_df.shape[0]//2 < int(thresh):
        print(f"{temp}/{fixed}_seed{seed}_ds{thresh}_chop.sam has less than {thresh} fragments")
    else:
        sam_df["QNAME"] = sam_df["temp"].apply(lambda row: row.split("\t")[0])  # take only read ID
        fragments = sam_df.drop_duplicates(subset=["QNAME"], inplace=False)  # remove duplicates/ convert reads to fragments (paired reads = 1 pair)
        fragments = fragments[["QNAME"]]
        # choose fragments to keep and subsample in pandas
        permuted = fragments["QNAME"].sample(n=int(thresh), random_state=int(seed))
        merged = sam_df.merge(permuted, how="right")
        merged = merged[["temp"]]
        merged.to_csv(f"{temp}/{fixed}_seed{seed}_ds{thresh}_tmp1.txt", index=False, header=False, sep="\t")
        # remove excess quotation marks from subsampled sam file
        sed_command = "sed 's/^.\{1\}//g' " + f"{temp}/{fixed}_seed{seed}_ds{thresh}_tmp1.txt | sed 's/.$//' > {temp}/{fixed}_seed{seed}_ds{thresh}_tmp2.txt"
        shell_command(sed_command, expand_vars=False)
        # put header back onto sam file
        shell_command(f"cat {temp}/{fixed}_seed{seed}_ds{thresh}_chop_header.txt {temp}/{fixed}_seed{seed}_ds{thresh}_tmp2.txt > {temp}/{fixed}_seed{seed}_ds{thresh}_final.sam")
        # convert sam to bam
        shell_command(f"samtools view -b -o {temp}/{fixed}_seed{seed}_ds{thresh}_final.bam {temp}/{fixed}_seed{seed}_ds{thresh}_final.sam")
        # copy to destination
        shell_command(f"rsync -avz {temp}/{fixed}_seed{seed}_ds{thresh}_final.bam {destination_jobID}/seed{seed}/downsample{thresh}/{fixed}.bam")
    # delete tmp files
    shell_command(f"rm {temp}/{fixed}_seed{seed}_ds{thresh}*")
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = elapsed_time//60, elapsed_time%60
    print(f"finished {destination_jobID}/seed{seed}/downsample{thresh}/{fixed}.bam; took {int(minutes)} minutes and {seconds:.2f} seconds \n", flush=True)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: <exe> <path to config file>', flush=True)
        print(len(sys.argv), sys.argv, flush=True)
        sys.exit(-1)
    else:
        _, config_path = sys.argv
        print(config_path, flush=True)
        with open(config_path, "r") as config:
            args = [a.split(":")[1].split("\n")[0] for a in config.readlines()]
            paths, fixnames_file, coverage_file, temp, destination, thresholds, seeds, jobID = args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7]
        print(paths, fixnames_file, coverage_file, temp, destination, thresholds, seeds, jobID, flush=True)

        # source, temp, and destination folders
        source = paths.split("/*")[0]
        if temp.endswith("/"): # remove end slashes if present
            temp = temp[:-1]
        if destination.endswith("/"):
            destination = destination[:-1]

        # fixnames file (determines input files)
        if fixnames_file == "":
            paths = [os.basename(p).split(".")[0] for p in glob.glob(paths)]
            fixnames_dict = dict(zip(paths, paths))
        else:
            fixnames_df = pd.read_csv(fixnames_file, sep="\t", header=None)
            fixnames_dict = dict(zip(fixnames_df[0].to_list(), fixnames_df[1].to_list()))
            paths = fixnames_df[0].to_list()

        # convert thresholds from string to list
        thresholds = thresholds.split(",")

        # convert seed(s) to list
        if len(seeds.split("-"))==2:
            seeds = [str(s) for s in range(int(seeds.split("-")[0]), int(seeds.split("-")[1]) + 1)]
        else:
            seeds = sorted(seeds.split(","))

        # create necessary folders
        for folder in [f"{destination}/{jobID}",f"{temp}/{jobID}"]:
            try:
                shell_command(f"mkdir {folder}")
            except:
                print(f"folder already exists: mkdir {folder}")
        temp = f"{temp}/{jobID}"

        for seed in seeds:
            try:
                shell_command(f"mkdir {destination}/{jobID}/seed{seed}")
            except:
                print(f"folder already exists: mkdir {destination}/{jobID}/seed{seed}")
            for thresh in thresholds:
                try:
                    shell_command(f"mkdir {destination}/{jobID}/seed{seed}/downsample{thresh}")
                except:
                    print(f"folder already exists: {destination}/{jobID}/seed{seed}/downsample{thresh}")

        # coverage file
        coverage_df = pd.read_csv(coverage_file, sep="\t")
        coverage_dict = dict(zip(coverage_df["sample"].to_list(), coverage_df["coverage"].to_list()))

        # generate job info
        all_thresh_jobs = sum([[[p,t] for t in thresholds] for p in paths],[])
        all_seed_jobs = tuple(sum([[(j[0],j[1],s) for s in seeds] for j in all_thresh_jobs], []))
        print(f"all jobs: {len(all_seed_jobs)}", flush=True)
        all_seed_jobs = [a for a in all_seed_jobs if os.path.exists(f"{destination}/{jobID}/seed{a[2]}/downsample{a[1]}/{fixnames_dict[a[0]]}.bam")==False]
        print(f"jobs to do: {len(all_seed_jobs)}", flush=True)

        # run jobs
        # PARALLELIZE IF DESIRED
        for f in all_seed_jobs:
            downsample(f[0], fixnames_dict[f[0]], coverage_dict[fixnames_dict[f[0]]], source, temp, f"{destination}/{jobID}", f[1], f[2])
        print("all files processed!", flush=True)




