# Detecting Meiotic Drive in Drosophila Virilis Pooled Sequence Data


This project attempts to re-implement the pooled sequencing methodology described in Wei et al. 2017 "[A Pooled Sequencing Approach Identifies a Candidate Meiotic Driver in Drosophila](https://pubmed.ncbi.nlm.nih.gov/28258181/)" 
on Drosophila Virilis pooled sequence data. 



## Usage

`af_calc.py` outputs a CSV file, per chromosome, containing the number of SNPs at biallelic segregating sites in the BC progeny that correspond to the v48 and v47 strains respectively. The outputted CSV files will be used by `plot_af.py` to produce meiotic drive diagnostic plots and output read depth plots.

```
python3 af_calc.py [-h] [-vcf VCF] [-vcf_index_file VCF_INDEX] 
                        [-bam BAM] [-bai BAI] [-o OUT_FILE]
```

Command-line arguments:
- `-vcf`: the VCF file containing all shared polymorphic sites between v48 and v47 
- `-vcf_index`: the index file for the VCF file
- `-bam`: the BAM file containing all the BC pooled sequences
- `-bai`: the index file for the BAM file
- `-o`: the name of the CSV file to which each chromosome's output will be written


`plot_af.py` creates meiotic drive diagnostic plots from the CSV files outputted from `af_calc.py`. It can also produce read depth plots that show the average read depth per site within a bin.
```
python3 plot_af.py [-h] [-csv CSV] [-contig CONTIG] [-window_size WINDOW_SIZE] 
                    [-fwer FWER] [-read_depth READ_DEPTH] [-spline SPLINE] 
                    [-o OUT_FILE]
```

Command-line arguments:
- `-csv`: name of the CSV file to perform the analysis on
- `-contig`: which chromosome / contig the CSV file contains data on
- `-window_size`: the size of the bins | default = 125,000
- `-fwer`: the family-wide error rate to use for the analysis | default = 0.01
- `-read_depth`: a boolean variable, if True read depth plots will be produced | default = False
- `-spline`: a boolean variable, if True a quintic spline will be fit to the SNP frequencies | default = False
- `-o`: the name of the outputted plot 


## Data Files

All data files can be found in this repository except the VCF and BAM files which can be found [here](https://drive.google.com/drive/folders/12b6tt0ZwQcZcno0uxiFK3C9PocLqxRFY?usp=sharing). 

I recommend just using `plot_af.py` on the provided CSV files, since `af_calc.py` takes an extremely long time to run.
Example commands to replicate the results provided in the report using `plot_af.py` are listed below:
```
python3 plot_af.py -csv Chr_2_pileup_allele_frequencies.csv -contig 2 -spline True -read_depth True 
python3 plot_af.py -csv Chr_3_pileup_allele_frequencies.csv -contig 3 -spline True -read_depth True 
python3 plot_af.py -csv Chr_4_pileup_allele_frequencies.csv -contig 4 -spline True -read_depth True 
python3 plot_af.py -csv Chr_5_pileup_allele_frequencies.csv -contig 5 -spline True -read_depth True 
python3 plot_af.py -csv Chr_6_pileup_allele_frequencies.csv -contig 6 -spline True -read_depth True 
python3 plot_af.py -csv Chr_X_pileup_allele_frequencies.csv -contig X -spline True -read_depth True
```

## Example Plot
![alt text](https://github.com/mohamedfaisa1/Detecting-Meiotic-Drive-in-Drosophila-Virilis-from-Pooled-Sequencing-Data/blob/main/data/Chr_4_allele_frequencies.png?raw=true)
