# RNA-SEQ

---

We describe steps involved in the analysis of RNASEQ.
This script has been developped for stranded pair-end RNASEQ.



---

1. Set up Environment
2. Set up ConfigFile for alignment only
3. Launch Alignment
4. Compute Statistics
5. Gene Expression Analysis


## Set up Environment

---

First you need to check that the following tools are installed.

_STAR_** : Aligner [here](https://github.com/alexdobin/STAR)

_Samtools_** : Bam handler [here](http://www.htslib.org/download/)

_wigToBigWig_** : Include in KentTools suite [here](http://hgdownload.soe.ucsc.edu/downloads.html#source_downloads)

Scripts available here are in Python 3. 
It's not required but advised to install Conda if python3 is not set up on your computer.
It will make things easier then for installing tools or switching to older python version if needed.

_Conda_** : [here](https://www.continuum.io/downloads)

## Set up ConfigFile.json for alignment only

---

You have a json file to config before doing your alignment.
The json file is a *key:value* listing which defines all parameters for the pipeline.
(number of cores used, path to output directory, path to input files in fastq.gz)


## Launch Alignment

---

It will launch an alignment with STAR.
It will creates a directory (path defined in your config file) with inside :
1. bam file
2. bigwig files per strand normalized by CPM (count per millions) for visualization purpose in UCSC.
3. Read count per gene - Sample_ReadsPerGene.tab (used for the gene expression step)
4. intermediate files used in the following steps of the pipeline.

```shell
python3 pathTo/align_V1.py -c pathToConfigFile/condition.json
```
You need to use this on each sample/replicate.

You have to manually grap each *read count per Gene* column to make only file.
Sort each file by geneName column before joining.
Depending on your library type , you should pick first or second column.
Select the output according to the strandedness of your data.
Here kit used for library preparation is of type directional first strand that means R2 is always forward.
The second read (read 2) is from the original RNA strand/template, first read (read 1) is from the opposite strand. 
Tcheck to get more details on the library preparation.[Here](https://www.youtube.com/watch?v=n2XEsw7EJLw&feature=youtu.be) 
[Here](http://chipster.csc.fi/manual/library-type-summary.html)

Select second column which name is *GoodPairs(F2R1+/F1R2-)*
At the end, it should like what follows :

| GeneName  | Condition1Replicat1    | Condition1Replicat2  | Condition2Replicat1  | ... |
| ---   |:---:| ---: | ----:| -----:|
| Gene1 | 100 | 140  | ...  |
| Gene2 | 0   |  0   | ...  |
| Gene3 | ...  | ... | ...  |

## Compute Statistics

---

This will create a file named *statistics* in you current directory with :

1. NameOfYourSample
2. Input Reads
3. Uniquely Mapped (%)
4. Multimapping (%)

```shell
python3 pathTo/computeStatistics.py -d pathToDirOfSample -n nameOfYourSample
```
You need to use this on each sample/replicate.


## Gene Expression Analysis

First, you need to create a csv file describing your experiments.
Column name can be changed but then you need to make changes directly inside code.
*Sample_id* and *condition* columns are mandatory.
It should looks as follows :

| sample_id  | condition    | group  |
| ---   |:---:| ---: | 
| Cond1rep1 | Cond1 | 2016  | 
| Cond1rep2 | Cond1   |  2016   | 
| Cond1rep1 | Cond1  |  2015  |

This will be your first input.
Then you will pass also the matrice containing all your reads count per gene and per sample.You did that the step before.
You have to give the names of the condition you want to compare between each others.

```R
diff_expRscript ${PATH_TO_SCRIPT}/diff_exp.R  --dir ${PATH_TO_DATA}/[DIR_NAME] --cond1 [COND1]  --cond2 [COND2]  ${PATH_TO_DATA}/[DESIGN.csv] ${PATH_TO_DATA}/[GENE_READ_COUNT.csv] > ${PATH_TO_DATA}/logR/[ANY_NAME].out
```
This script will give you list of differentially expressed genes with fold changes and associated pvalues. Graphics ( PCA plot, Heatmap, Histogram , ...) that describes quality of your data will also be created.
One file with all pvalue and fc will be created per comparison.

You can merge all this files into one. You need first to sort them by gene ensembl id column to be sure the fusion will be correct.

Note : There is also others text files created which already apply filter on the foldchange(log2FC > 1.5) .You can change that directly in the code if you want.

```R
merged_all_FC_files ${PATH_TO_SCRIPT}/merged_all_FC_files.R  --dir ${PATH_TO_DATA}/[DIR_OUTPUT]  
```
You need to modify inside the code the path to the files because they are hard-coded.
At the end, it will create a file (*FC_collapse_with_pval.csv*) in the directory you passed as parameter.
