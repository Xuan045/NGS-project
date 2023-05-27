# NGS-project

## Variant calling pipelines

![image](https://github.com/Xuan045/NGS-project/assets/86905456/f3237df6-f467-4eba-9ffc-53208bae551f)

1. GATK best practices
  - Data preprocessing was performed using GATK and Picard tools. For more detailed information, please refer to [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery)
  - SNVs and indels were called using GATK HaplotypeCaller in GVCF mode, followed by joint analysis. For more information, please refer to [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) 

2. [Google DeepVariant](https://github.com/google/deepvariant)

  - Per sample variants were called in GVCF mode.
  - Per sample GVCF files were merged using [GLnexus](https://github.com/dnanexus-rnd/GLnexus).
  - We follow [Best practices for multi-sample variant calling with DeepVariant (WES trio demonstration)](https://github.com/google/deepvariant/blob/r0.9/docs/trio-merge-case-study.md) to perform the multi-sample variant calling.
  - Reference: Yun, T., Li, H., Chang, P., Lin, M. F., Carroll, A., & McLean, C. Y. (2021). Accurate, scalable cohort variant calls using DeepVariant and GLnexus. Bioinformatics, 36(24), 5582-5589. https://doi.org/10.1093/bioinformatics/btaa1081

## Annotation
- We used ANNOVAR for variant annotation.
- For more information about ANNOVAR, please refer to the [official document](https://annovar.openbioinformatics.org/en/latest/)
- Reference: Wang, K., Li, M., & Hakonarson, H. (2010). ANNOVAR: Functional annotation of genetic variants from high-throughput sequencing data. Nucleic Acids Research, 38(16), e164. https://doi.org/10.1093/nar/gkq603
