# GRENE-net Phase 1

This repository contains the code and analysis used in the paper "Rapid adaptation and extinction across climates in synchronized outdoor evolution experiments of *Arabidopsis thaliana*" (bioRxiv: https://www.biorxiv.org/content/10.1101/2025.05.28.654549v1).

## Abstract 

Climate change is threatening species with extinction, and rapid evolutionary adaptation may be their only option for population rescue over short ecological timescales. However, direct observations of rapid genetic adaptation and population dynamics across climates are rare across species. To fill this gap, we conducted a replicated, globally synchronized evolution experiment with the plant *Arabidopsis thaliana* for 5 years in over 30 outdoor experimental gardens with distinct climates across Europe, the Levant, and North America. We performed whole-genome sequencing on ∼70,000 surviving reproductive individuals and directly observed rapid and repeatable adaptation across climates. Allele frequency changes over time were parallel in experimental evolution replicates within the same climates, while they diverged across contrasting climates—with some allele frequency shifts best explained by strong selection between −46% to +60%. Screening the genome for signals of rapid climate adaptation identified a polygenic architecture with both known and novel adaptive genetic variants connected to important ecological phenotypes including environmental stress responses, CAM5 and HEAT SHOCK FACTORs, and germination and spring flowering timing, CYTOCHROME P450s and TSF. We found evolutionary adaptation trends were often predictable, but variable across environments. In warm climates, high evolutionary predictability was associated with population survival up to 5 years, while erratic trends were an early warning for population extinction. Together, these results show rapid climate adaptation may be possible, but understanding its limits across species will be key for biodiversity forecasting.


## Repository Structure

The repository is organized into several key directories:

```
.
├── Genome_Environment_Association/     # Contains all GEA-related analyses
│   ├── binomial_regression/           # GLM with binomial link for GEA
│   ├── climate_distance/             # Climate distance calculations
│   ├── follow_up_sign_blocks/        # Follow-up analysis of significant blocks
│   ├── gwas/                        # Climate GWAS studies
│   ├── kendall_tau/                 # Kendall's tau correlation analyses for GEA
│   ├── lfmm/                        # Latent factor mixed model analyses for GEA
│   ├── manhattan_plot_GEA_annotated/ # Manhattan plot generation for GEA
│   ├── preproc_scripts/             # Data preprocessing scripts
│   ├── signficant_intersection_GEA_models/ # Intersection of significant GEA results
│   └── wza/                         # Weighted Z-score analysis for P-value aggregation across haplotypes
├── Local_adaptation/                 # Analysis of local adaptation patterns
├── Population_Evolution/             # Overall population evolution analyses
├── Predictability/                   # Analysis of predictability of adaptation
├── demography_analysis/              # Population demographic analyses
├── genomic_offset/                   # Genomic offset calculations
├── genomic_patterns/                 # Analysis of genomic wide patterns over space and time
└── temporal_dynamics/                # Analysis of temporal patterns
```

## Data Availability

The input data files required to run these analyses will be made available through Zenodo soon 


