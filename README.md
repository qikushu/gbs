## Overview

Genotyping-by-sequencing (GBS) is a protocol used for the discovery and genotyping of SNPs. 
ddRAD-seq involves massively analyzing the terminal sequences of genomic DNA fragments that have been cut with restriction enzymes, and by comparing these sequences between varieties and strains, DNA polymorphisms (differences in DNA sequences) can be identified.

This repository provides easy-to-use scripts that change VCF files produced by Tassel-GBS pipeline v2 (Glaubitz et al. 2014) into ABH genotype format. The ABH format produced is ready for use in QTL analyses with the R/qtl package (). It also includes extra scripts that help running the LB-impute software to fill in missing genotypes in Hidden Markov Model approach (Fragoso et al. 2016).

## Prerequisites

Before you begin, ensure you have met the following requirements:

- VCF files outputted by the Tassel-GBS pipeline v2.
- Parents of a biparental cross included in the VCF files.
- LB-impute is available.
- Tassel pipepline is available.

## Instruction

### vcf to hapmap (Without imputation)
VCF file `in.vcf` is converted to a hapmap format file `out.hmp.txt` by the Tassel pipeline.

```bash
$TASSELPIPELINE="$HOME/tassel-5-standalone/run_pipeline.pl"
perl $TASSELPIPELINE -Xmx5g -fork1 -vcf in.vcf -export out -exportType Hapmap -runfork1
```

The resulting `out.hmp.txt` is converted to the ABH genotype format. The script `hmp_to_rqtl_cp.pl` is included in this repository.

```bash
less out.hmp.txt | perl hmp_to_rqtl_cp.pl out.hmp.txt Parent_For_A Parent_For_B OUT
```

The first and second argument values of the perl script are the parental variety 
names in the hapmap file (Parent_For_A Parent_For_B). 
The third argumet value directs to include (IN) or not to include (OUT) parental genotypes to the R/QTL format file.
The fourth and consequent argument value are the comma-delimited values of individuals excluded from output.
converstion from hapmap to r/qtl format.

This script produces the following four files. 
```
joinmap_CP.txt
joinmap_NNNP.txt
joinmap_LMLL.txt
joinmap_HKHK.txt
ABH.txt
```

### vcf to hapmap (With imputation)
Missing genotypes in a VCF file `in.vcf` is imputed by LB-impute software.

```bash
less in.vcf | perl prepareLBimpute.pl > new.vcf

java -Xmx2g -jar ../toolbox/LB-Impute/LB-Impute.jar -method impute -offspringimpute -f new.vcf -recombdist 10000000 -window 5 -o imputed.vcf -parents Parent_A, Parent_B
```

The resulting VCF file `imputed.vcf` is converted to hapmap and ABH genotype file.


## References
Bradbury, P.J., Z. Zhang, D.E. Kroon, T.M. Casstevens, Y. Ramdoss, and E.S. Buckler (2007) TASSEL: software for association mapping of complex traits in diverse samples. Bioinformatics 23: 2633–2635.

Broman, K. W., Wu, H., Sen, S., & Churchill, G. A. (2003). R/qtl: QTL mapping in experimental crosses. Bioinformatics (Oxford, England), 19(7), 889–890.

Fragoso, C. A., Heffelfinger, C., Zhao, H., & Dellaporta, S. L. (2016). Imputing Genotypes in Biallelic Populations from Low-Coverage Sequence Data. Genetics, 202(2), 487–495. 

Glaubitz, J.C., T.M. Casstevens, F. Lu, J. Harriman. R.J. Elshire, Q. Sun, E.S. Buckler (2014) TASSEL-GBS: A High Capacity Genotyping by Sequencing Analysis Pipeline. PLoS ONE 9(2): e90346. 


## Publication

