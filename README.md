# S4 Polygenic Risk Score Pipeline

This repository provides the **S4 programs** and scripts for generating and evaluating polygenic risk scores (PRS).  
It contains the C++ source code, executables, and example workflows, along with pre-configured runs and documentation.

The tutorial below walks through the step-by-step workflow for preparing references, selecting SNPs, shrinking models, generating PRS, and running downstream analyses.  


## Repository Structure

The repository includes the following folders:

- **BootstrapTesting/**  
  Scripts for estimating the standard deviation of accuracy differences between polygenic models.

- **FinnGenRuns/**  
  Example execution runs for model testing using summary statistics from the **FinnGen** study.

- **LDpredRuns/**  
  Scripts and configurations for generating **LDpred** models.

- **PGMfiles/**  
  A collection of the best-performing polygenic models, optimised using either Biobank or FinnGen datasets.

- **S4_runs/**  
  Sample runs showcasing the usage of the **S4 model**.

- **Testing_Method/**  
  Documentation describing the model testing procedures.

- **doc/**  
  General documentation for the project.

- **exe/**  
  Linux executable files (mostly statically compiled, with the exception of dynamic linking).

- **inc/**  
  Header files used in the C++ source code.

- **lib/**  
  Libraries required by the C++ program.

- **src/**  
  Source code files written in **C++**.

---

## Workflow

This tutorial documents the step-by-step workflow for running the **S4 programs** to generate and evaluate polygenic risk scores (PRS).

The pipeline covers reference preparation, SNP selection, model shrinking, PGS generation, and downstream statistical analysis.



### **Step 0. Prerequisites**

* [BGENIX](https://enkre.net/cgi-bin/code/bgen/doc/trunk/doc/wiki/bgenix.md)
* [QCTOOL](https://www.chg.ox.ac.uk/~gav/qctool_v2/documentation/download.html)
* `gawk`
* GNU `parallel` (optional, depends on the HPC system in your institute)


### **Step 1. Prepare Reference Files**

Generate reference genotype files in **BGEN format** for your reference samples. The `haplotypes_chr*.bgen` files are derived from phased variant data from the 1000 Genomes Project (1KGP) high-coverage autosomal data. Download the files from the following FTP sites:


* GRCh38 (hg38): [1000G High-Coverage Data](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/)

* GRCh37 (hg19): [1000G Phase 3 Data](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)

To download VCF and index files for chromosomes 1 to 22 (GRCh38), run:

```bash
for i in {1..22}; do 
  wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
  wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
done

```

**File Conversion**

1. **Split into Haplotypes**: Extract haplotypes from the downloaded VCF files.

2. **Convert to .gen Format**: Use `awk` to convert haplotypes into genotype files (`.gen`).

3. **Convert to BGEN Format**: Use `qctool` to convert `.gen` files to `haplotypes_chr*.bgen` format for downstream analysis.


**Note:**
* Index Files: The `.tbi` index files are required for efficient processing in downstream steps.
* Storage: Ensure sufficient disk space, as high-coverage VCF and BGEN files are large.

### **Step 2. Create BGENIX Index Files**

Index BGEN files for fast access:

```bash
parallel --dryrun bgenix -g referencefile_chr{}.bgen -index ::: {1..22} > BgenixRuns.sh
```

### **Step 3. Create SNP List Files**

List SNPs for one chromosome:

```bash
bgenix -g haplotypes_chr1.bgen -list > haplotypes_snps_chr1.txt
```

List SNPs for multiple chromosomes in parallel:

```bash
parallel --dryrun 'bgenix -list -g referencefile_chr{}.bgen > reference_snps_chr{}.txt' ::: {1..22} > Create_snp_list_files.sh
```

### **Step 4. Prepare Summary Statistics**

Prepare per-chromosome summary statistics with the following **tab-separated columns**:

1. SNP Name
2. Position
3. Effect Allele
4. Reference Allele
5. Frequency
6. R²
7. Odds Ratio
8. Standard Error
9. P-value
10. Functional score (optional, if using functional annotations)

```bash
1:100:C:A	100	A	C	0.183418	0.625168	0.0181689	0.0190609	0.0206728 3.212
1:101:G:A	101	A	G	0.513657	0.685011	0.0385628	0.0208292	0.0398586 1.681
1:102:A:G	102	G	A	0.283657	0.690986	-0.0131923	0.0197965	0.367652 2.351
1:103:A:T	103	T	A	0.196798	0.770082	-0.0129992	0.0197811	0.230602 8.296
1:104:T:G	104	G	T	0.526066	0.680632	-0.0136653	0.0197622	0.353208 6.203
```

### **Step 5. Calculate N-value**

Approximate the effective sample size (**N-value**) using chromosome 2 (largest #variants):

```bash
gawk '$8>0 && $5>0.1 && $5<0.9 && $6>0.8{a=1/($5*(1-$5)*$8^2);print a}' summary_statistics_chr2.txt | sort -k1gr | head -10000 | tail -1
```

* Filter:

  * `se > 0`
  * `0.1 < MAF < 0.9`
  * `imputation R² > 0.8`
* Round the resulting value to **two significant figures**.
* Example: `27689.2 → 28000`

### **Step 6. Run S4\_Select**

Run per chromosome to select SNPs:

```bash
S4_select --functional --ancestry --effect --correct --corrdamp 0.99 \
-m 0.005 -r 0.05 --maxr 0.85 --maxindr 0.95 --dist 3 --r2impute 0.3 -z 5 \
--excludesnps ambindel --genmap genetic_map_hg38_chr1.txt \
-d haplotypes_chr1.bgen --refsnpfile shared_b38_snps_chr1.txt \
--reg 0.05 -a b38_ancall_phenotype_chr1.gz b38_eur_phenotype_chr1.gz b38_eas_phenotype_chr1.gz b38_afr_phenotype_chr1.gz \
-i select_haplo_continental.txt -c 1 -p 1e-4 -n 120000 \
> phenotype_noambindel_r85_1e-4_b38_chr1_snps.txt
```

**Key parameters:**

* `--functional`: Use functional annotations (optional).
* `--ancestry`: Enable multi-ancestry mode (optional).
* `--effect`: Use z-score changes for grouping.
* `--correct`: Correct for imputation quality.
* `--corrdamp`: Damping factor for correlations (default: 0.99).
* `-m 0.005`: MAF threshold.
* `-r 0.05`: R² threshold for grouping (default: 0.05).
* `--maxr 0.85`: Max correlation (corrected).
* `--maxindr 0.95`: Max individual correlation.
* `--dist 3`: Max genetic distance.
* `--r2impute 0.3`: Minimum imputation accuracy.
* `-z 5`: QC threshold.
* `--excludesnps ambindel`: Excludes ambiguous and indel SNPs.
* `--genmap`: Specifies the genetic map file.
* `-d`: Reference file for calculating LD.
* `--refsnpfile`: List of SNPs to consider.
* `--reg 0.05`: Regularization value to dampen estimates (default: 0.05).
* `-a`: Ancestry files.
* `-i`: File containing individuals or ancestral groups to include (optional).
* `-c`: Chromosome number.
* `-p 1e-4`: P-value/R² threshold for SNP inclusion.
* `-n`: N-value from Step 5.

**Note:** 
* Omit `--ancestry` and `-i` flags if your GWAS summary statistics include only one ancestry group.
* For multiple ancestry groups, list GWAS summary statistics files in sequence, as shown above.

### **Step 7. Concatenate SNP Files**

Merge all chromosome outputs:

```bash
cat phenotype_noambindel_r85_1e-4_b38_chr{1..22}_snps.txt > phenotype_noambindel_r85_1e-4_b38_snps.txt
```

### **Step 8. Run S4\_Shrink**

Shrink SNPs into a polygenic model:

```bash
S4_shrink --functional -w select_haplo_continental.txt --corrdamp 0.99 --reg 0.05 --approx 0.01 \
-i phenotype_noambindel_r85_1e-4_b38_functional_snps.txt -n 10000 -r haplotypes_chr \
--correct --nvalue 120000 --seed 1245 -a 0.1 -b 0.7 -p 5e-6 \
-o s4_phenotype_anc_r85_1e-4_functional_0.1_0.7_5e-6_pgm.gz
```

**Key parameters:**

* `--functional`: Use functional annotations (optional).
* `-w`: Enable multi-ancestry mode.
* `--corrdamp`: Damping factor for correlations (default: 0.99).
* `--reg 0.05`: Regularization value (default: 0.05).
* `--approx 0.01`: Speeds up calculations by skipping high-penalty SNPs.
* `-i`: SNP file generated by S4_select.
* `-n 10000`: Number of MCMC permutation runs.
* `-r`: Prefix for reference LD files.
* `--correct`: Corrects for imputation quality.
* `--nvalue`: N-value calculated in Step 5.
* `--seed 1245`: Sets a random seed for reproducibility.
* `-a 0.1 -b 0.7 -p 5e-6`: Parameters for the MCMC model.
* `-o`: Output file containing the polygenic model.

### **Step 9. Generate Polygenic Scores**

```bash
PGS_generate -d haplotypes_chr{1..22}.bgen \
-w s4_phenotype_anc_r85_1e-4_functional_0.1_0.7_5e-6_pgm.gz \
-o s4_phenotype_anc_r85_1e-4_functional_0.1_0.7_5e-6_weightings.txt
```

**Key parameters:**

* `-d`: List og BGEN files for PGS calculation.
* `-w`: PGM file generated by S4_shrink.
* `-o`: Output file containing PGS scores.

### **Step 10. Analyse PGS Results**

```bash
MultiPGSstats -g CC --pheno example_phenotypes.txt \
-w s4_phenotype_anc_r85_1e-4_functional_0.1_0.7_5e-6_weightings.txt
```

**Key parameters:**

* `-g CC`: Phenotype to analyse.
* `--pheno`: Phenotype file.
* `-w`: Weighting file generated by PGS_generate

**Output columns:**

* Unadjusted OR
* Standard error of raw OR
* Likelihood ratio test chi²
* AUC
* SE of raw values
* OR per SD
* 95% CI lower bound
* 95% CI upper bound

### **Step 11. Summarize Statistics**

```bash
MultiPGSsumstats -s b38_finngen_r12_breast.gz \
-b haplotypes_chr0.bgen \
-r shared_b38_snps_chr0.txt \
-w s4_phenotype_anc_r85_1e-4_functional_0.1_0.7_5e-6_pgm.gz \
--genmapfile genetic_map_hg38_chr0.txt --dist 3 \
-i select_haplo_eur.txt \
> s4_phenotype_anc_r85_1e-4_functional_0.1_0.7_5e-6_finngen_stats.txt
```

**Key parameters:**

* `-s`: Summary statistics file.
* `-b`: Haplotype reference files.
* `-w`: PGM weightings.
* `--genmapfile`: Genetic map files.
* `--dist 3`: Genetic distance threshold for SNP correlations.
* `-i`: File specifying included individuals based on the ethnicity of samples.

### **Step 12. Store Correlations (Optional)**

For larger PGMs, use CorrMat to store correlations and speed up calculations.

```bash
CorrMat -b haplotypes_chr1.bgen -r shared_b38_snps_chr1.txt \
-w s4_phenotype_anc_r85_1e-4_functional_0.1_0.7_5e-6_pgm.gz \
--genmapfile genetic_map_hg38_chr1.txt --dist 3 \
-i select_haplo_eur.txt \
-o eur_corr_phenotype_anc_example_chr1.gz
```

Use stored correlations in `MultiPGSsumstats`:

```bash
MultiPGSsumstats -s b38_biobank_r12_pheno.gz \
-c eur_corr_pheno_anc_example_chr0.gz \
-w s4_phenotype_anc_r85_1e-4_functional_0.1_0.7_5e-6_pgm.gz \
> s4_phenotype_anc_r85_1e-4_functional_0.1_0.7_5e-6_precalculated_corr_biobank_stats.txt
```

---

## Notes

* [Genetic map files](https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/) (`--genmap`) are required (hg19 or hg38). Please download it and separate by chromosomes.
```bash
$ head genetic_map_hg38_chr1.txt
chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)
1 55550 0.0000000000 0.000000000000000
1 632942 0.0000000000 0.000000000000000
1 633147 0.0000000000 0.000000000000000
1 785910 2.6858076690 0.410292036939447
1 788439 2.8222713027 0.417429561063975
1 788511 2.9813105581 0.417644215424158
1 792862 2.9806151254 0.430612871834774
1 794568 3.0780969498 0.435864105231133
1 805477 3.0751332930 0.469410734324470
```
* `-i` include files specify individuals/ancestral groups. Omitting this can produce empty outputs.
* Recommended to test multiple parameter sets (`-p`, `-a`, `-b`, `phi`) and compare results using AUC or chi².

---

## Citation

If you use this repository, please cite the following works:

1. Dareng, E.O., Tyrer, J.P., Barnes, D.R. *et al.*  
   Polygenic risk modeling for prediction of epithelial ovarian cancer risk.  
   *European Journal of Human Genetics* 30, 349–362 (2022).  
   [https://doi.org/10.1038/s41431-021-00987-7](https://doi.org/10.1038/s41431-021-00987-7)

2. Tyrer, J.P., Peng, P.-C., DeVries, A.A. *et al.*  
   Improving on polygenic scores across complex traits using select and shrink with summary statistics (S4) and LDpred2.  
   *BMC Genomics* 25, 878 (2024).  
   [https://doi.org/10.1186/s12864-024-10706-3](https://doi.org/10.1186/s12864-024-10706-3)

3. Baierl, J., Tyrer, J.P., Lai, P.-H. *et al.*  
   S4-Multi: enhancing polygenic score prediction in ancestrally diverse populations.  
   *medRxiv* 2025.01.24.25321098 (2025).  
   [https://doi.org/10.1101/2025.01.24.25321098](https://doi.org/10.1101/2025.01.24.25321098)
