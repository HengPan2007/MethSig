Introduction to the MethSig pacakge
================
Heng Pan
2020-07-26

## Setup

``` r
setRepositories(ind = 1:2)
install.packages("devtools")
devtools::install_github("HengPan2007/MethSig", build_vignettes = T)
library(MethSig)
```

## Introduction

MethSig is used to infer DNA methylation drivers from bisulfite
converted data of cancer cohorts.

## Prerequisite

### 1\. Data alignment and processing

Bisulfite converted data is aligned and processed as described in the
published book chapter (Pan *et al.*, *Cancer Systems Biology*, 2018).

### 2\. Covariates matrix preparation

A matrix of covariates used in the beta regression model is needed. The
matrix needs to include a column of hugo gene symbol and at least one
column of a covariate. An example can be loaded by *invisible(matCV)*.

``` r
matCV <- invisible(matCV)
```

This example matrix includes Hugo symbol (Hugo), average promoter DHcR
level in control samples (DHcR\_Normal), average promoter PDR level in
control samples (PDR\_Normal), average gene expression level in control
samples (GEXP\_Normal) and DNA replication time (Reptime). Users can
define their own covariates.

    #>     Hugo DHcR_Normal  PDR_Normal GEXP_Normal    Reptime
    #> 1 A4GALT   -0.126499  0.12442306  1.33315372 -1.1379561
    #> 2   AAAS   -0.126499 -0.67038461  0.97661449 -1.1462431
    #> 3   AACS   -0.126499 -0.37628615  0.64288254 -0.5661545
    #> 4  AADAT   -0.126499 -0.08258285 -0.04805863  1.6091775
    #> 5  AAGAB   -0.126499 -0.61130311  2.51903184 -0.5081457
    #> 6   AAK1   -0.126499 -0.28766984  0.64981415 -0.5537241

### 3\. Files containing information of hypermethylated cytosines (HCs)

A tab-separated values input file (without a header line) contains
details of differentially methylated cytosines (DMCs) with following
columns (V1 to V11): chr, pos, numC in control, numC + numT in control,
numC in tumor, numC + numT in tumor, CpG methylation ratio (tumor
methylation / control methylation), chi-squared test p-value, adjusted
p-value, significance, hyper or hypo in tumor. Details of generating
this type of files were described in Pan *et al.*, *Cancer Systems
Biology*, 2018. An example file can be found in
extdata/DMC.SRR2069925.txt. Notably, the name needs to be in
DMC.sample\_name.txt format.

    #>     V1     V2  V3  V4 V5 V6    V7      V8 V9 V10 V11
    #> 1 chr1  10497 461 562 79 90 1.071 0.23325  1   -   -
    #> 2 chr1  10525 521 563 83 90 0.997 1.00000  1   -   -
    #> 3 chr1 136876 600 635 34 36 1.001 1.00000  1   -   -
    #> 4 chr1 136895 546 639 32 36 1.044 0.74231  1   -   -
    #> 5 chr1 713376  41  44 16 19 0.911 0.51846  1   -   -
    #> 6 chr1 713388  38  46 18 19 1.145 0.37194  1   -   -

### 4\. Files containing information of proportion of discordant reads (PDR) at CpGs

A tab-separated values input file contains details of single CpG PDR
with following columns: chr, start, strand, ConMethReadCount,
ConUMethReadCount, DisReadCount, NAReadCount. Details of PDR were
described in Landau *et al.*, *Cancer Cell*, 2014.
extdata/pdrCall\_from\_Bismark.py can be used to call PDR of single CpG
from Bismark (Krueger *et al.*, *Bioinformatics*, 2011) output files
starting with CpG\_OB or CpG\_OT. An example file can be found in
extdata/PDR.SRR2069925.txt. Notably, the name needs to be in
PDR.sample\_name.txt
    format.

    #>    chr  start strand ConMethReadCount ConUMethReadCount DisReadCount
    #> 1 chr1  10497      +                0                 0            0
    #> 2 chr1  10525      +                0                 0            0
    #> 3 chr1 136876      +                0                 0            0
    #> 4 chr1 136895      +                0                 0            0
    #> 5 chr1 713376      +                0                 0            0
    #> 6 chr1 713388      +                0                 0            0
    #>   NAReadCount
    #> 1          90
    #> 2          90
    #> 3          36
    #> 4          36
    #> 5          19
    #> 6          19

## Input matrix generation

### 1\. Promoter DHcR calculation

Promoter (defined as ± 2kb window centered on RefSeq transcription start
site) hypermethylation is measured using differentially hypermethylated
cytosine ratio (DHcR), defined as the ratio of hypermethylated cytosines
(HCs) to the total number of promoter CpGs profiled. HCs of each sample
are defined as CpGs at which DNAme is statistically higher than the
average DNAme of control samples (false discovery rate = 20%,
Chi-squared test). Only CpGs with read depth greater than 10 are
included in the analysis. RRBS data of normal samples are used as
control. An implemented function *makeHG19Promoters* can be used to
provide promoter annotations of hg19 RefSeq genes. Users and define
their own
annotation.

``` r
dhcr <- promoterDHcR(file_name = system.file("extdata", "DMC.SRR2069925.txt", package = "MethSig"),
                     pro = makeHG19Promoters())
head(dhcr)
```

    #>       Hugo    Depth HCs CpGs  DHcR
    #> 41  ABCB10 44.91818   0  110 0.000
    #> 62   ABCD3 38.69118   0   68 0.000
    #> 93    ABL2 36.35714   0   28 0.000
    #> 113  ACADM 32.36364   0   33 0.000
    #> 120  ACAP3 40.36500   5  200 0.025
    #> 123  ACBD3 36.91892   0   37 0.000

### 2\. Promoter PDR calculation

If all the CpGs on a specific read are methylated, or all the CpGs on a
read are unmethylated, the read is classified as concordant; otherwise
it is classified as discordant. At each CpG, the PDR is equal to the
number of discordant reads that cover that location divided by the total
number of reads that cover that location. The PDR of promoter is given
by averaging the values of individual CpGs, as calculated for all CpGs
within the promoter of interest with equal or greater than 10 reads
covering at least 4
CpGs.

``` r
pdr <- promoterPDR(file_name = system.file("extdata", "PDR.SRR2069925.txt", package = "MethSig"),
                   pro = makeHG19Promoters())
head(pdr)
```

    #>       Hugo        PDR
    #> 41  ABCB10 0.06439566
    #> 58   ABCC8 0.33333333
    #> 62   ABCD3 0.02900665
    #> 88    ABI1 0.04192547
    #> 93    ABL2 0.02704556
    #> 103  ABTB2 0.04380624

### 3\. Input matrix generation

As mentioned in **Prerequisite** section, users need to put
DMC.sample\_name.txt and PDR.sample\_name.txt files in the input\_dir
folder. Also, a user defined covarites matrix is needed.

``` r
ds <- makeInputMatrix(names_list = as.list("SRR2069925"),
                matCV = invisible(matCV),
                pro = makeHG19Promoters(),
                input_dir = system.file("extdata", "", package = "MethSig"))
head(ds)
```

    #>           Id   Hugo DHcR_Normal  PDR_Normal GEXP_Normal    Reptime  PDR_Tumor
    #> 1 SRR2069926 A4GALT   -0.126499  0.12442306  1.33315372 -1.1379561 0.11196453
    #> 2 SRR2069926   AAAS   -0.126499 -0.67038461  0.97661449 -1.1462431 0.10410410
    #> 3 SRR2069926   AACS   -0.126499 -0.37628615  0.64288254 -0.5661545 0.11056835
    #> 4 SRR2069926  AADAT   -0.126499 -0.08258285 -0.04805863  1.6091775 0.13897421
    #> 5 SRR2069926  AAGAB   -0.126499 -0.61130311  2.51903184 -0.5081457 0.05119315
    #> 6 SRR2069926   AAK1   -0.126499 -0.28766984  0.64981415 -0.5537241 0.06419811
    #>   Depth_Tumor CpGs_Tumor DHcR_Tumor
    #> 1    58.05263         57 0.01754386
    #> 2    24.27273         11 0.00000000
    #> 3    32.40541         37 0.00000000
    #> 4    40.22917         48 0.02083333
    #> 5    30.52174         23 0.08695652
    #> 6    46.37500         48 0.02083333

### Sample-specific hypermethylation inference

Expected promoter DHcR of tumor samples is estimated by beta regression
model and expected DHcR is tested against observed DHcR to infer
hypermethylation
status.

``` r
pval <- pvalueBetaReg(formula = as.formula("DHcR_Tumor_Beta~DHcR_Normal+PDR_Normal+GEXP_Normal+Reptime+PDR_Tumor+Depth_Tumor+CpGs_Tumor"),
                      data = invisible(inputMat))
head(pval)
```

    #>           Id   Hugo DHcR_Normal  PDR_Normal GEXP_Normal    Reptime  PDR_Tumor
    #> 1 SRR2069926 A4GALT   -0.126499  0.12442306  1.33315372 -1.1379561 0.11196453
    #> 2 SRR2069926   AAAS   -0.126499 -0.67038461  0.97661449 -1.1462431 0.10410410
    #> 3 SRR2069926   AACS   -0.126499 -0.37628615  0.64288254 -0.5661545 0.11056835
    #> 4 SRR2069926  AADAT   -0.126499 -0.08258285 -0.04805863  1.6091775 0.13897421
    #> 5 SRR2069926  AAGAB   -0.126499 -0.61130311  2.51903184 -0.5081457 0.05119315
    #> 6 SRR2069926   AAK1   -0.126499 -0.28766984  0.64981415 -0.5537241 0.06419811
    #>   Depth_Tumor CpGs_Tumor DHcR_Tumor DHcR_Tumor_Beta Beta_Response
    #> 1    58.05263         57 0.01754386    1.758860e-02    0.07016196
    #> 2    24.27273         11 0.00000000    4.636499e-05    0.06188360
    #> 3    32.40541         37 0.00000000    4.636499e-05    0.07261282
    #> 4    40.22917         48 0.02083333    2.087777e-02    0.10148911
    #> 5    30.52174         23 0.08695652    8.699482e-02    0.04220170
    #> 6    46.37500         48 0.02083333    2.087777e-02    0.05943839
    #>   Beta_Precision    Pvalue
    #> 1       3.079141 0.4468126
    #> 2       3.079141 0.8068584
    #> 3       3.079141 0.8557167
    #> 4       3.079141 0.5635852
    #> 5       3.079141 0.1427462
    #> 6       3.079141 0.3724881

### Cohort-prevalent hypermethylation inference

Wilkinson p-value combination method is used to determine if promoter
hypermethylation is over-represented in the cohort. To eliminate the
effect of cohort size on p-value combination results, MethSig randomly
samples equal number of patients iteratively and uses lower quartile of
combined p-values to infer hypermethylation.

``` r
pval <- pvalueCombine(data = invisible(pvalByGenePt))
head(pval)
```

    #>              Hugo Rank Sample_Size       Pvalue      Padjust
    #> HIST1H3E HIST1H3E    1          23 2.109371e-28 2.300269e-24
    #> PCDHGB7   PCDHGB7    2          18 7.069669e-28 3.210187e-24
    #> CARTPT     CARTPT    3          22 8.831325e-28 3.210187e-24
    #> TMPRSS12 TMPRSS12    4          22 3.647620e-27 8.078755e-24
    #> ELF3         ELF3    5          24 3.704152e-27 8.078755e-24
    #> CYP4F8     CYP4F8    6          12 5.184899e-27 8.224423e-24
