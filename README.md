# metagenome_matcher
Verify sample identities and detect processing errors using host/human SNPs inferred from metagenomic sequencing. Use one or both of the following approaches depending on your dataset: 
1) Metagenome v. Genotype: compare SNPs inferred from metagenomic sequencing of a sample to independently-obtained SNPs (such as SNPs from a genotype microarray) to identify the true donor of a sample. Requires a separate source of SNP data in addition to metagenomic sequencing.
2) Metagenome v. Metagenome: compare SNPs inferred from metagenomic sequencing between samples. Requires multiple metagenomic samples per individual. 

## IMPORTANT: File Naming Conventions
The following naming conventions are critical for determining which donors supplied which sample(s) according to sample labels.

All metagenomic samples must be labeled in the following convention: donorname_samplename. For example, sample A1 from study participant Bob could be labeled Bob_A1. If samples are not already in this format, they can be renamed using the "fixnames" function of preprocess.py (see below). If supplying a PLINK fileset for Approach #1, set both FID and IID as the donorname. For example, Bob would appear in the .fam file of the PLINK fileset with FID=Bob and IID=Bob. 

## Notes on Config Files
```preprocess.py```, ```metagenome_v_genotype.py```, ```metagenome_v_metagenome.py```, and ```downsample.py``` all intake config ```.txt``` files where the user specifies parameters. 
- on each line of config files, type parameters directly after the colon (do not place a space after the colon).
- when supplying full paths to directories, exlude the final slash (e.g, ```/my/full/path``` instead of ```/my/full/path/```)
- you will supply paths to both output folders and temp folders for intermediate files that will be deleted. These folders do not have to exist (the program will create them if they do not already exist), but the folder they are to be created in must exist.
- choose job IDs that are unique and descriptive, as they will be used to label both temp and final output files and folders
- supply any necessary tables or dataframes headerless, index-free, space-delimited text files

## Prerequisites and Installation:
### Operating System Requirements
Metagenome Matcher was designed for use on macOS and Linux, and tested on the following:
- macOS: Ventura 13.2.1 
- Linux: Red Hat Enterprise Linux 9.3 (Plow)

### Python Packages (versions tested in parentheses)
Tested on Python 3.11.8, 3.13.5
- numpy (1.24.4, 2.3.2)
- pandas (2.3.1, 2.3.2)
- matplotlib.pyplot (3.7.2, 3.10.0)
- seaborn (0.12.2, 0.13.2)

### Other Programs
- samtools (1.17, 1.22.1)
- bcftools (1.21, 1.22)
- PLINK (1.90b6.26, 1.9.0-b.7.7)

### Installation
Install in seconds with: ```git clone https://github.com/korem-lab/metagenome_matcher.git```

## Preprocessing
Regardless of your chosen approach, all samples will need processing with preprocess.py. Input files must be in ```.bam``` or ```.bam.gz``` format (alignments of metagenomic sequencing data to the host/human genome). All preprocessing parameters are specified in the config file (see ```preprocess_config.txt```). Output (processed samples in samtools mpileup format, as well as a a ```.tsv``` file with coverage of the host genome for each sample) are stored in the specified destination folder.

### Notes on parameters
- To rename  files in the final output, supply a space-delimited file with original file names in the first column and new names in the second column.
- If you plan to use Approach #1, ensure that the chromosome names in your .bam files are the same as the chromosome names your additional genotype dataset. If they differ (e.g., "CM000663.2" vs "1"), you can supply a space-delimited file: the first column of the file should be chromosome names as they appear in sequencing files, and the second column should be chromosome names as they appear in your SNP array genotypes. If chromosome names match between sequencing and SNP array genotypes, or if you will be comparing sequencing to sequencing (mpileup swapped analysis), leave this blank.
- You must supply a list of loci to compare SNPs (space-separated file with chromosome names in first column, position coordinates in second column). If you plan to use Appraoch #1, ensure that chromosome names in this list match chromosome names in SNP array genotypes. Comparing only at select loci can greatly reduce the size of each file and thus speed up downstream comparisons and decrease storage constraints. If using Approach #1, a good choice would be to filter for only loci in profiled by the genotype microarray. ```target_loci.txt``` within this repo, a selection of the loci profiled by the Infinium Multi-Ethnic Global Array (Illumina, USA), can be used.

### On the Command Line: 
```python preprocess.py preprocess_config.txt```

## Approach 1: Metagenome v. Genotype Comparisons
Compare host SNPs inferred from metagenomic sequencing to genotypes from an independent source (e.g., microarray genotyping). metagenome_v_genotype.py requires preprocessed mpileup files (generated via preprocess.py) as input for metagenomic sequencing data. For the independently-obtained SNPs, please supply the path to PLINK fileset (e.g., ```/my/full/path/files``` for ```files.bed```, ```files.bim```, and ```files.fam``` within the directory ```/my/full/path/files```). All parameters, including the paths to input files, are specified in the config file (see ```metagenome_v_genotype_config.txt```). Output (a ```.tsv``` with the number of SNPs with data in both the metagenome and array file and the percentage of SNPs concurrent between data sources, as well as graphs, if requested in the config file) is stored in the specified destination folder. An individual is likely a sample's donor if metagenome vs. genotype comparison yields a lax concurrence score of 98% or more.

### Modes:
- "all": compare each sample to all individuals with available genotype data.
- "self": compare each sample to the genotypes of its labeled donor (confirm high concurrence between labeled sample donor and sample)
- [path to pairs file]: supply a space-delimited file with specific metagenome:genotype pairs to compare (each pair on its own line)

### On the Command Line:
python ```metagenome_v_genotype.py metagenome_v_genotype_config.txt```

## Approach 2: Metagenome v. Metagenome Comparisons
Compare host SNPs inferred from metagenomic sequencing between samples. metagenome_v_metagenome.py requires preprocessed mpileup files (generated via ```preprocess.py```) as input. All parameters, including the path to input files, are specified in the config file (see ```metagenome_v_metagenome_config.txt```). Output (a ```.tsv``` with the number of SNPs with sequencing coverage in both metagenomes and the percentage of SNPs concurrent between metagenomes, as well as graphs, if requested in the config file) is stored in the specified destination folder. Two samples are likely from the same donor if metagenome vs. metagenome comparison yields a strict concurrence score of 90% or more.

### Modes:
- "all": compare all possible pairs of samples within the specified input directory
- "same": compare only samples labeled as being from the same donor (and thus expected to be similar)
- [path to a pairs file]: supply a space-delimited file with specific pairs of samples to compare (each pair on its own line)
- [path to .fam file]: only samples with labeled donors NOT in the specified .fam file (samples whose labeled donors lack independently-obtained genotypes and are thus ineligible for analysis with Approach #1) will be analyzed. These samples will be compared with other samples from their labeled donor.

### On the Command Line:
```python metagenome_v_metagenome.py metagenome_v_metagenome_config.txt```

## Downsampling
To evaluate the performance of metagenome_matcher when faced with low sequencing coverage, we downsampled our files to specific readcounts (of the human genome) using downsample.py. At each readcount threshold, we conducted 10 experiments. Our downsampling script requires prior knowledge of the original samples' coverage of the host genome (supplied as a ```.tsv``` with each sample followed by its coverage on its own line) which can be generated with ```preprocess.py```.

As with preprocessing, if you desire to rename your files in the final output, you can supply a .tsv with original and new file names (one original:new name pair per line). The number of experiments conducted per readcount threshold is specified via the number of random seeds given in the config file  (see ```downsample_config.txt```). All other parameters (including readcount thresholds) are also specified in the config file. Output downsampled files are stored in the specified destination folder. 

### On the Command Line:
```python downsample.py downsample_config.txt```

## Demo
Simulated data (both ```.bam``` files and array genotypes in PLINK format) are provided for a demo. The Metagenome Matcher python commands (steps 5, 6, and 7) should each complete within one second.
1. Download the human reference Genbank assembly (GCA_000001405.29) as a fasta from: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
2. In the same directory as the .fna file that you downloaded ```GCA_000001405.29_GRCh38.p14_genomic.fna```, run: ```samtools faidx GCA_000001405.29_GRCh38.p14_genomic.fna```
3. Navigate to the cloned ```metagenome_matcher``` directory.
4. Modify the following ```.txt``` config files with your own paths: ```toy_preprocess_config.txt```, ```toy_metagenome_v_metagenome_config.txt```, and ```toy_metagenome_v_genotype_config.txt```
5. Run ```python preprocess.py toy_preprocess_config.txt```
6. Run ```python metagenome_v_metagenome.py toy_metagenome_v_metagenome_config.txt```. Compare the results output to ```metagenome_matcher/toy_data/toy_results/SNPsleuth_toymvm``` to those in ```metagenome_matcher/toy_data/toy_expected_results/SNPsleuth_toymvm```
7. Run ```python metagenome_v_genotype.py toy_metagenome_v_genotype_config.txt```. Compare the results output to ```metagenome_matcher/toy_data/toy_results/SNPsleuth_toymvg``` to those in ```metagenome_matcher/toy_data/toy_expected_results/SNPsleuth_toymvg```
