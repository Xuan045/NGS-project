# NGS-project

## Variant calling pipelines

![image](https://github.com/Xuan045/NGS-project/assets/86905456/8ea4fb56-4aa5-4021-94f0-91e0ada75c75)

1. GATK best practices
  - Data preprocessing was performed using GATK and Picard tools. For more detailed information, please refer to [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)
  - SNVs and indels were called using GATK HaplotypeCaller in GVCF mode, followed by joint analysis. For more information, please refer to [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) 

2. [Google DeepVariant](https://github.com/google/deepvariant)

  - Per sample variants were called in GVCF mode.
  - Per sample GVCF files were merged using Picard MergeVcfs.
  


