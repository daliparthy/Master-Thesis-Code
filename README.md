# Solid-State Dewetting Microstructural & Crystallographic Analysis

Automated Python and MATLAB scripts developed for quantitative Scanning Electron Microscopy (SEM) image segmentation and Electron Backscatter Diffraction (EBSD) crystallographic analysis. This codebase supports the Master's thesis investigating the influence of annealing atmospheres on the solid-state dewetting (SSD) of platinum (Pt) thin films on sapphire substrates across multiple film thicknesses at the **Chair for Micro- and Nanostructure Research, Friedrich-Alexander-Universität Erlangen-Nürnberg (FAU)**.

---

## Overview of Workflows

### 1. SEM Image Segmentation & Particle Analysis (Python)
Batch-processes top-view SEM micrographs to quantify dewetting evolution from continuous/perforated films to isolated islands:
* **Automated Thresholding:** Applies Adaptive/Otsu/Triangle thresholding algorithms to separate dewetted Pt islands from the exposed sapphire substrate.
* **Morphological & Geometric Metrics:** Computes covered area fraction, equivalent circular diameter, particle circularity.
* **Batch Reporting:** Aggregates particle-level and image-level statistics across experimental conditions and exports structured `.xlsx` workbooks for statistical fitting and plotting in OriginPro.

### 2. EBSD Grain Reconstruction & Texture Analysis (MATLAB / MTEX)
Batch-processes EBSD orientation maps (`.ang` files) to evaluate grain growth, crystallographic texture, and grain boundary character:
* **Data Filtering & Grain Reconstruction:** Imports `.ang` datasets, filters low-quality indexing points by Confidence Index (CI), and reconstructs grains for face-centered cubic (FCC) platinum.
* **Grain Size Statistics:** Calculates number-weighted and area-weighted grain size distributions and pools datasets across multiple scans for representative statistics.
* **Texture & Boundary Characterization:** Computes Orientation Distribution Functions (ODFs), quantifies texture strength via the Texture Index ($J$-index), and analyzes Coincidence Site Lattice (CSL) grain boundary distributions.

---

## Packages & Dependencies

### Python Environment
Tested with **Python 3.9 or later**. Required third-party packages:

| Package | Import Name | Purpose in Pipeline |
| :--- | :--- | :--- |
| **OpenCV** (`opencv-python`) | `cv2` | Image preprocessing, Otsu/Triangle thresholding, contour detection, and morphological filtering |
| **NumPy** | `numpy` | Array operations, equivalent circular diameter ($\text{ECD}$), and volume calculations |
| **pandas** | `pandas` | Tabular data structuring, multi-sample pooling, and statistical aggregation |
| **openpyxl** | `openpyxl` | Exporting multi-sheet particle and summary statistics to Excel (`.xlsx`) |

Standard library modules used: `os`, `pathlib`, `glob`, `math`.

### MATLAB Environment
Tested with **MATLAB R2026a**. Required toolbox:

| Toolbox | Purpose in Pipeline |
| :--- | :--- |
| **MTEX** (version 6.1.1) | `.ang` file import, CI filtering, FCC Pt crystal symmetry definition, grain reconstruction, ODF calculation, J-index quantification, and CSL boundary analysis |

---

## Installation & Setup

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/daliparthy/Master-Thesis-Code.git](https://github.com/daliparthy/Master-Thesis-Code.git)
   cd Master-Thesis-Code
