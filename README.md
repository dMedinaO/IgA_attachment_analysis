# Immune system up-regulation in securely attached infants during early childhood

This repository contains the source files and supplementary information for the implementations and use cases presented in the work:

Jorge González-Puelma<sup>1,2†</sup>, Lindybeth Sarmiento Varón<sup>2†</sup>, Jessica Vidal<sup>3</sup>, Constanza Ceroni<sup>2</sup>, Sebastián Escobedo<sup>2,4</sup>, Roberto Uribe-Paredes<sup>4,5</sup>, David Medina-Ortiz<sup>4,5</sup>, Rodrigo A. Cárcamo<sup>3*</sup>, Marcelo A. Navarrete<sup>1,2*</sup><br>

Immune system up-regulation in securely attached infants duringearly childhood. <br>https://doi.org/XXXXXX<br>

<sup>*1*</sup><sub>Escuela de Medicina, Universidad de Magallanes, Avenida Bulnes 01855, Punta Arenas, Chile.</sub> <br>
<sup>*2*</sup><sub>Centro Asistencial Docente e Investigación, Universidad de Magallanes, Av. Los Flamencos 01364, Punta Arenas, Chile.</sub> <br>
<sup>*3*</sup><sub>Facultad de Psicología y Humanidades, Universidad San Sebastián, Gral Lagos 1163, 5110693 Valdivia, Chile.</sub> <br>
<sup>*4*</sup><sub>Departamento de Ingeniería en Computación, Universidad de Magallanes, Avenida Bulnes 01855, Punta Arenas, Chile.</sub> <br>
<sup>*5*</sup><sub>Centre for Biotechnology and Bioengineering, CeBiB, Universidad de Chile, Beauchef 851, Santiago, Chile.</sub> <br>
<sup>*†*</sup><sub>These authors contributed equally </sub> <br>
<sup>*\**</sup><sub>Corresponding author</sub> <br>

---
## Table of Contents
- [A summary of the proposed work](#summary)
- [Folder description](#folder)
- [Requirements](#requirements)
- [Raw and processed data](#raw-data-and-preprocessing)
- [General comments](#comments)
- [Contacts](#contacts)
---

<a name="summary"></a>

## Immune system up-regulation in securely attached infants during early childhood

Early caregiving relationships shape the coordination of stress and immune systems, yet their biological correlates in early infancy remain insufficiently understood. This study examined whether attachment relationships are associated with mucosal immune function and stress physiology. Thirty-five infants (mean age = 16.6 months) were classified as securely or insecurely attached using the Strange Situation Procedure. Salivary secretory IgA (SIgA) was collected in the morning and afternoon at two time points, and cumulative cortisol was quantified from hair samples.

Securely attached infants showed higher morning SIgA concentrations and more stable intra-day immune profiles compared with insecurely attached children. No group-level differences were observed for cumulative cortisol, but immune–endocrine associations revealed that higher cortisol was linked to lower morning SIgA and greater intra-day fluctuation. Bayesian regression models supported consistent directional effects, and machine-learning analyses confirmed that SIgA-based features accurately predict attachment type.

Our findings support the idea that secure attachment fosters stable coordination between immune and endocrine systems during a critical stage of early development. These effects could be readily captured at very early stages of life, identifying SIgA as a potential biomarker of early socioemotional environments. By integrating behavioral, immunological, and computational approaches, this study provides evidence for the biological embedding of attachment and highlights the potential of non-invasive biomarkers to support early identification of psychosocial vulnerability.

---

<a name="requirements"></a>

## Requirements

This project was implemented using Python 3.12. In the file [requirements.txt](requirements.txt) you will find all modules and dependencies required to run the experiments.

Based on the requirements, an environment file [environment.yml](environment.yml) was generated. To crete the environment, please write the following command line in the linux terminal:

```
conda env create -f environment.yml
```

With the environment created, you can activate the environment

```
conda activate iga_attachment
```

<a name="raw-data-and-preprocessing"></a>

## Raw data and preprocessing

The raw data of this paper is available in the folder: [raw_data](raw_data/).

The following content is available:

- [cortisol](raw_data/cortisol/raw_data_measures.csv): Data measures for cortisol
- [sIgA](raw_data/mortality_data.xlsx): Data for mortality analysis

The folders contain the following information:

- [data_for_figures](data_for_figures/): Data necessary to generate the figures integrated in the paper.
- [figure_generation](figure_generation/): Notebooks employed to generate the figures.
- [figures_for_paper](figures_for_paper/): Figures generated for the paper
- [generated_mmodels](generated_mmodels/): Models generated during the work
- [inference_process](inference_process/): Notebooks implemented to generate the inference models and the predictive methods
- [preprocessing_data](preprocessing_data/): Notebooks implemented to process the raw data.
- [processed_data](processed_data/): Processed datasets from the raw data
- [raw_data](raw_data/): Raw data employed in this work
- [results_process](results_process/): Statistical results using the processed datasets
- [statistical_analysis](statistical_analysis/): Notebooks implemented to assess statistical analysis.
- [tmp_figures](tmp_figures/): Temporal figures generated for subpanels.

<a name="contacts">Contacts</a>

- Corresponding author: Marcelo Navarrete [marcelo.navarrete@umag.cl](marcelo.navarrete@umag.cl)
- Corresponding author: Rodrigo Cárcamo [rodrigo.carcamol@uss.cl](rodrigo.carcamol@uss.cl)
- First author: Jorge González-Puelma [jorge.gonzalez@umag.cl](jorge.gonzalez@umag.cl)
- First author: Lindybeth Sarmiento-Varón [lindybeth.sarmiento@umag.cl](lindybeth.sarmiento@umag.cl)
- Data anlysis/implementation: David Medina-Ortiz [david.medina@umag.cl](david.medina@umag.cl )

## Citation  

If you use this repository, please cite:  

> González-Puelma, J., Sarmiento Varón, L., Vidal, J., Ceroni, C., Escobedo, S., Uribe-Paredes, R., Medina-Ortiz, D., Cárcamo, R. A., & Navarrete, M. A. (2025). mmune system up-regulation in securely attached infants during early childhood. [Journal name], [volume(issue)], [pages]. https://doi.org/[DOI].  
