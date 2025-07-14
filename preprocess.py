from Utils import shell_command
import os
import re
import glob
import sys
import pandas as pd
import time
import random

def filter_mp(sample, tmp, chrom_dict, merge_col, ji):
    chrom_dict = dict(zip(chrom_dict[0], chrom_dict[1]))
    totals = pd.DataFrame(columns=["chromosome", "position", "reference", "coverage", "readbases", "rqual", "mqual",
                                   "Merge"])
    ACTG = {'A', 'C', 'T', 'G'}
    docu = f"{tmp}/{ji}_{sample}.pileup.gz"  # path to mpileup
    headers = ["chromosome", "position", "reference", "coverage", "readbases", "rqual", "mqual"]  # mpileup headers
    print("begin loading mpileup", docu, flush=True)
    with pd.read_csv(docu, names=headers, sep="\t", quoting=3, chunksize=(10 ** 6)) as reader:
        for mpile in reader:
            if mpile.empty:
                print(f"{sample} empty", flush=True)
            else:
                mpile = mpile[mpile["reference"].apply(lambda ref: ref.upper() in ACTG)]  # get rid of entries without reference
                mpile["chromosome"] = mpile["chromosome"].apply(
                    lambda x: chrom_dict[x] if x in chrom_dict.keys() else "SKIP")  # chromosome decoding
                mpile["Merge"] = mpile["chromosome"] + "_" + mpile["position"].astype(str)  # Merge label for mpileup
                merged = mpile.merge(merge_col, how="inner", on="Merge")
                totals = pd.concat([totals, merged])
    final = merge_col.merge(totals, how="inner", on="Merge")
    savepath = f"{tmp}/{ji}_{sample}_final.pileup.gz"
    final.to_csv(savepath, sep="\t", compression="gzip")

def pipeline(sample, sample_fixedname, source, tmp, dest, chrom_dict, merge_col, regf, ref, rlen, gz_opt, ji):
    start_time = time.time()
    print(f"starting preprocessing for {sample}", flush=True)

    if os.path.exists(f"{dest}/mpileup_filtered/{sample_fixedname}.pileup.gz"):
        print(f"{sample} mpileup already done", flush=True)
        return

    # STEP1: copy files over and sort alignments, decompressing first if necessary
    if gz_opt == 1:
        file_stats = os.stat(f"{source}{sample}.bam.gz")
        if file_stats.st_size == 0:
            print(f"{sample} file is empty", flush=True)
            return
        copy_cmd = f"rsync -avz {source}{sample}.bam.gz {tmp}/{ji}_{sample}.bam.gz"
        print(copy_cmd, flush=True)
        shell_command(copy_cmd)
        gunzip_cmd = f"gunzip -c {tmp}/{ji}_{sample}.bam.gz > {tmp}/{ji}_{sample_fixedname}.bam"
        print(gunzip_cmd, flush=True)
        shell_command(gunzip_cmd)
        shell_command(f"rm -f {tmp}/{ji}_{sample}.bam.gz") # remove COMPRESSED version
        print(f"removed compressed file: rm -f {tmp}/{ji}_{sample}.bam.gz", flush=True)
    else:
        file_stats = os.stat(f"{source}{sample}.bam")
        if file_stats.st_size == 0:
            print(f"{sample} file is empty", flush=True)
            return
        copy_cmd = f"rsync -avz {source}{sample}.bam {tmp}/{ji}_{sample_fixedname}.bam"
        print(copy_cmd, flush=True)
        shell_command(copy_cmd)
    sample = sample_fixedname # from now on, we use the fixed name
    sort_cmd = f"samtools sort -o {tmp}/{ji}_{sample}_sorted.bam {tmp}/{ji}_{sample}.bam"
    print(sort_cmd, flush=True)
    shell_command(sort_cmd)

    # STEP2: index sorted alignments
    index_cmd = f"samtools index {tmp}/{ji}_{sample}_sorted.bam"  # index goes into same folder
    print(index_cmd, flush=True)
    shell_command(index_cmd)

    # STEP2.5: calculate and record coverage
    covg_calc_cmd1 = f"samtools depth {tmp}/{ji}_{sample}_sorted.bam | "
    covg_calc_cmd2 = "awk '{sum+=$3} END { print sum/" + str(rlen) + "}' > " # assumes specific reference genome
    covg_calc_cmd3 = f"{tmp}/coverage/{ji}_{sample}.txt"
    covg_calc_cmd = covg_calc_cmd1 + covg_calc_cmd2 + covg_calc_cmd3
    if os.path.exists(f"{tmp}/coverage/{ji}_{sample}.txt")==False: #file does not exist
        shell_command(covg_calc_cmd, expand_vars=False)
        print(covg_calc_cmd, flush=True)
    elif os.path.getsize(f"{tmp}/coverage/{ji}_{sample}.txt")==0: #file exists but is empty
        shell_command(covg_calc_cmd, expand_vars=False)
        print(covg_calc_cmd, flush=True)
    else:
        print(f"coverage for {sample} already calculated",flush=True)

    # STEP3: filter indexed, sorted alignment file to include only SNPs in the array
    filter_cmd = f"samtools view -bh -L {regf} {tmp}/{ji}_{sample}_sorted.bam > {tmp}/{ji}_{sample}_filtered.bam"
    print(filter_cmd, flush=True)
    shell_command(filter_cmd)

    # STEP4: sort filtered alignment
    sort_filtered_cmd = f"samtools sort -o {tmp}/{ji}_{sample}_filtered_sorted.bam {tmp}/{ji}_{sample}_filtered.bam"
    print(sort_filtered_cmd, flush=True)
    shell_command(sort_filtered_cmd)

    # STEP5: index sorted, filtered alignment
    index_filtered_cmd = f"samtools index {tmp}/{ji}_{sample}_filtered_sorted.bam"
    print(index_filtered_cmd, flush=True)
    shell_command(index_filtered_cmd)

    # STEP6: generate pileup
    mpileup_cmd = f"samtools mpileup -f {ref} -s {tmp}/{ji}_{sample}_filtered_sorted.bam | gzip -9 -cf > {tmp}/{ji}_{sample}.pileup.gz"
    print(mpileup_cmd, flush=True)
    shell_command(mpileup_cmd)

    # STEP7: extra filtering on locus level
    print("ready to filter: ", sample, tmp, flush=True)
    filter_mp(sample, tmp, chrom_dict, merge_col, ji)

    # STEP 8: copy final file over
    copy_mp_cmd = f"rsync -avz {tmp}/{ji}_{sample}_final.pileup.gz {dest}/mpileup_filtered/{sample}.pileup.gz"
    print(copy_mp_cmd)
    shell_command(copy_mp_cmd)

    # STEP 9: remove all remaining tmp files for specific sample
    rm_bam_cmd = f"rm {tmp}/{ji}_{sample}*"
    print(rm_bam_cmd, flush=True)
    shell_command(rm_bam_cmd)

    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes, seconds = elapsed_time // 60, elapsed_time % 60
    print(f"finished {sample}: took {int(minutes)} minutes and {seconds:.2f} seconds \n", flush=True)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: <exe> <path to config file>', flush=True)
        print(len(sys.argv), sys.argv, flush=True)
        sys.exit(-1)
    else:
        _, config_path = sys.argv
        print(config_path, flush=True)
        with open(config_path,"r") as config:
            args = [a.split(":")[1].split("\n")[0] for a in config.readlines()]
            paths, fixnames, tmp_path, dest_inp, chrom_path, targets, reference, ref_len, job_id = args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8]
        print(paths, fixnames, tmp_path, dest_inp, chrom_path, targets, reference, ref_len, job_id, flush=True)
        pathnoslash = [i[:-1] if i.endswith("/") else i for i in [tmp_path, dest_inp]] # remove end slash if included
        tmp_path, dest_inp = pathnoslash[0], pathnoslash[1]

        if paths[-2:] == "gz": # if .bams are compressed (.bam.gz)
            gz_option = 1
        else: # .bams are uncompressed
            gz_option = 0

        source_path = paths.split("*")[0] # path to root folder with input files, everything before name of file
        dest_path = f"{dest_inp}/preprocess_{job_id}" # name of folder within destination folder that will contain all output
        rei = "^" + source_path + "(.*)" + paths.split("*")[1] + "$"
        basenames = [re.match(fr"{rei}", x).group(1) for x in glob.glob(paths)] # get base names of all input files (exludes extensions and path to folder)
        
        if fixnames!="": # fixnames file provided, change sample names as desired
            fixnames = pd.read_csv(fixnames, sep="\s+", header=None, names=["original", "fixed"]) # read in .tsv
            fixnames = fixnames.set_index('original')['fixed'].to_dict() # convert to dict; og names = key, new names = vals
        else: # no fixnames file provided, use sample names as-is
            fixnames = dict(zip(basenames, basenames))
        # only include samples within given fixnames file (all if no file provided) AND not already completed (check if output exists)
        basenames = [b for b in basenames if b in fixnames and os.path.exists(f"{dest_path}/mpileup_filtered/{fixnames[b]}.pileup.gz")==False]
        print(source_path, dest_path, basenames, len(basenames), flush=True)

        #make folders necessary. if folders already exist, will keep them
        to_mkdir = [dest_inp, dest_path, f"{dest_path}/mpileup_filtered", tmp_path, f"{tmp_path}/coverage"]  # figure out which folders you need
        for m in to_mkdir:  # make each folder
            try:
                shell_command(f"mkdir {m}")
                print(f"{m}: folder created", flush=True)
            except:
                print(f"{m}: folder already exists", flush=True)

        # take list of desired loci and create region file and merge column
        targets_source = "/".join(targets.split("/")[:-1])
        targets = pd.read_csv(targets, sep="\s+", header=None, dtype=str)
        targets["Merge"] = targets[0] + "_" + targets[1]
        merge_column = targets[["Merge"]]
        if chrom_path == "":
            chromosome_dict = (tuple(targets[0].unique()), tuple(targets[0].unique()))
        else:
            chrom_df = pd.read_csv(chrom_path, sep="\s+", dtype=str, header=None)  # read in chromosome name conversion dictionary
            chromosome_dict = (tuple(chrom_df[0].to_list()), tuple(chrom_df[1].to_list()))
        rev_dict = dict(zip(chrom_df[1].to_list(), chrom_df[0].to_list()))
        targets["seq_chrom"] = targets[0].apply(lambda row: rev_dict[row])
        targets["start"] = targets[1].apply(lambda row: int(row) - 1 if int(row) > 1 else row)
        targets["end"] = targets[1].apply(lambda row: int(row) + 1)
        reg_file = targets[["seq_chrom", "start", "end"]]
        regf_path = f"{targets_source}/region_file_{job_id}.bed"
        reg_file.to_csv(regf_path, header=False, index=False, sep=" ")

        for b in basenames: # if desired, modify code to parallelize this step. Otherwise, each file will be processed sequentially.
            pipeline(b, fixnames[b], source_path, tmp_path, dest_path, chromosome_dict, merge_column, regf_path, reference, ref_len, gz_option, job_id)
        print("finished preprocessing all files", flush=True)

        covg_files = []
        for i in glob.glob(f"{tmp_path}/coverage/{job_id}_*.txt"):
            covg_df = pd.read_csv(i, sep="\t", header=None, names =["coverage"])
            covg_df["sample"] = i.split(f"/{job_id}_")[-1].split(".")[0]
            covg_files.append(covg_df)
        covg_cat = pd.concat(covg_files)
        covg_cat = covg_cat[["sample", "coverage"]]  # reorder columns
        covg_cat.to_csv(f"{tmp_path}/coverage/coverage.tsv", sep="\t", index=False)
        shell_command(f"rsync -avz {tmp_path}/coverage/coverage.tsv {dest_path}/")
        shell_command(f"rm -rf {tmp_path}/coverage")
        print("finished creating coverage file", flush=True)
