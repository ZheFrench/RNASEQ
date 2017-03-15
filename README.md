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

## Set up ConfigFile for alignment only

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
Depending on your library type , you should pick first or second column.
Select the output according to the strandedness of your data.
Here kit used for library preparation is of type directional first strand that means R2 is always forward.
The second read (read 2) is from the original RNA strand/template, first read (read 1) is from the opposite strand. 
Tcheck to get more details on the library preparation.[here](https://www.youtube.com/watch?v=n2XEsw7EJLw&feature=youtu.be) 
Here[http://chipster.csc.fi/manual/library-type-summary.html)

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
1.nameOfYourSample
2.Input Reads
3.Uniquely Mapped (%)
4.Multimapping (%)

```shell
python3 pathTo/computeStatistics.py -d pathToDirOfSample -n nameOfYourSample
```
You need to use this on each sample/replicate.


## Gene Expression Analysis

