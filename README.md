
**Environmental filtering and anthropogenic pressures structure reef-water bacterial microbiomes across Red Sea coral reefs**

This repository contains the data and R scripts used in the manuscript:

**Environmental filtering and anthropogenic pressures structure reef-water bacterial microbiomes across Red Sea coral reefs**

Submitted to **Microbiome** (under revision).

**Authors**

Gloria Lisbet Gil Ramos, Eva Aylagas, João Cúrdia, Warren R. Francis, Karla Gonzalez, Micaela Sofia Dos Santos Justo, Andreia S. F. Farinha, Carolina Bocanegra-Castano, Diego Lozano-Cortés, Michael Berumen, and Susana Carvalho.


**Abstract**: 

Background: Reef-water microbiomes are increasingly recognized as integrative indicators of coral reef ecosystem condition because they rapidly respond to environmental change while reflecting biological, chemical, and physical processes occurring across the reef. However, the ecological drivers governing reef-water bacterial microbiome assembly across interacting natural and anthropogenic gradients remain poorly resolved, particularly at regional spatial scales. Here, we characterized reef-water bacterial microbiomes collected immediately above sixteen coral reefs spanning three regions of the central Saudi Arabian Red Sea. By integrating microbial community profiles with comprehensive environmental, benthic, and anthropogenic datasets, we investigated how chemical, physical, biological, and human-related drivers shape bacterial diversity, taxonomic composition, predicted functional profiles, and the abundance of ecologically important taxa.

Results: Reef-water bacterial microbiomes exhibited pronounced regional taxonomic turnover despite comparatively conserved predicted functional profiles, suggesting substantial functional redundancy across reefs. Multivariate analyses, Random Forest models, and Generalized Additive Models consistently identified sea surface temperature, distance to shore, chlorophyll-a, benthic habitat characteristics, nutrient availability, trace metals, and local anthropogenic pressure as the principal determinants of microbiome composition and diversity. Dominant bacterial taxa, including Prochlorococcus marinus and members of the order Rhodobacterales, displayed distinct responses to environmental gradients, revealing taxon-specific ecological niches and sensitivities to multiple stressors. Together, these findings demonstrate that reef-water microbiomes integrate environmental and anthropogenic signals across spatial scales and respond predictably to ecosystem condition.

Conclusions: Our findings show that reef-water bacterial microbiomes are structured through the combined effects of environmental filtering and local anthropogenic pressures, resulting in strong taxonomic turnover while maintaining comparatively stable predicted functional profiles. By identifying the principal drivers of reef-water microbiome assembly across one of the world's largest coral reef systems, this study advances our understanding of coral reef microbial ecology and highlights the potential of reef-water microbiomes as scalable, high-resolution indicators of ecosystem condition. These results provide a robust framework for incorporating microbiome-based approaches into coral reef monitoring, conservation, and ecosystem management under accelerating environmental change


**Repository contents**
This repository contains the data and scripts required to reproduce the analyses presented in the manuscript.

**Data**

The data are provided as .csv files and include:

Metadata associated with each water sample
Rarefied OTU tables used for diversity analyses and community composition.
Non-rarefied OTU tables used for Random Forest and Generalized Additive Model (GAM) analyses.
Model_GAMSandRF.csv, containing the predictor variables used to fit the Random Forest and GAM models.

**Scripts**

The repository includes the following R scripts:

Random_Forest.R – Random Forest analyses used to identify the main predictors of bacterial diversity, community composition, and selected bacterial taxa.
GAMS.R – Generalized Additive Models (GAMs) used to investigate nonlinear relationships between environmental variables and bacterial taxa/diversity metrics.

**Reproducibility**

All analyses were performed in R. (vversion 4.3.3). The scripts are intended to reproduce the statistical analyses presented in the manuscript using the data provided in this repository.

**Citation**

If you use the data or scripts in this repository, please cite:

Gil Ramos GL et al.
Environmental filtering and anthropogenic pressures structure reef-water bacterial microbiomes across Red Sea coral reefs.
Microbiome. (under revision)

