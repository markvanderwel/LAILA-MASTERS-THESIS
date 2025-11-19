# Master’s Thesis – Diversity Effects on Biomass Accumulation

This repository contains the reproducible analytical workflow for my master’s thesis, which investigates how taxonomic, functional, and phylogenetic diversity influence biomass accumulation in Atlantic Forest restoration sites.

The repository includes raw and processed spreadsheets and R scripts. Large raw data, GIS files, and heavy intermediates are intentionally excluded through the `.gitignore`.

---

## 📁 Repository Structure

### **01 Datasets/**
This folder contains the datasets used in the analyses.

01 Datasets/
├── 01_raw_data/ # Raw field measurements (species, traits, biomass, etc.)
├── 02_processed_data/ # Cleaned and standardized datasets used in the models
└── 04_original_data/ # Original file from fieldwork measurements, not to work with


---

### **02 R Scripts/**
All R scripts used in the thesis analyses.

02 R Scripts/
├── 01_phylo_master_analysis.R # Full phylogenetic diversity workflow
└── 02_functional_master_analysis.R # Full functional diversity workflow

(...) 


#### Master Scripts Description

- **01_phylo_master_analysis.R**  
  Complete pipeline for phylogenetic analyses  
  (tree construction, PD/MPD/MNTD, SES metrics, decoupled dissimilarities, PCPS, and more).

- **02_functional_master_analysis.R**  
  Complete functional diversity workflow  
  (trait processing, standardization, CWM, FDis, PCA, decoupled effects, diversity–productivity models).

Each master script consolidates the final, validated components of the analyses. Whenever an analysis is updated, the final code is added to these master files.There are some other files with intermediate steps.

---

## 📄 Project Files

- **README.md** — Project documentation  
- **MASTERS-THESIS.Rproj** — RStudio project file  
- **.gitignore** — Excludes heavy or sensitive files  

---

## 🔒 Data Availability

Large files (rasters, climate layers, GIS shapefiles, and heavy raw datasets) are excluded and therefore not tracked in this repository. These can be provided upon request.

---

## 📌 Notes

This repository is designed to:

- Support reproducibility of all statistical analyses  
- Maintain transparent workflow organization  
- Allow advisors and collaborators to access code and processed data  

---

