
```markdown
# SyNQA: Synergistic Network QUBO Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

## About SyNQA
Synergistic Network QUBO Analysis (SyNQA) is a computational framework that reframes biomarker discovery from a conventional univariate ranking to an interaction-aware combinatorial optimization problem. 

By formulating disease-specific interaction rewiring as an energy minimization problem of a Quadratic Unconstrained Binary Optimization (QUBO) model, SyNQA effectively extracts robust, non-redundant, and mechanistically interpretable feature subsets without relying on black-box machine learning algorithms. While originally applied to human gut microbiome data for colorectal cancer, its data-agnostic nature allows for potential adaptation to various high-dimensional omics modalities.

## Repository Structure
```text
.
├── data/
│   ├── demonstration/          # Dummy/demonstration datasets for testing
│   │   └── generate_demo_data.py # Script to generate dummy datasets
│   └── raw/                    # (Empty) Directory for user's full datasets
├── results/                    # Output directory for results and figures
├── SyNQA.py                    # Main execution script
├── environment.yml             # Mamba/Conda environment configuration
├── requirements.txt            # Python dependencies (pip)
├── LICENSE                     # MIT License
└── README.md
```

## Requirements and Installation
The pipeline is implemented in Python 3.10+. The optimization process utilizes Simulated Annealing provided by D-Wave's `dwave-neal` package. We strongly recommend using mamba (or conda) to set up the environment.

**1. Clone this repository:**
```bash
git clone [https://github.com/your-username/SyNQA.git](https://github.com/your-username/SyNQA.git)
cd SyNQA
```

**2. Create and activate the environment:**
```bash
mamba env create -f environment.yml
mamba activate synqa_env
```
> **Note:** If you prefer using pip, you can run `pip install -r requirements.txt` instead.

## Usage

### 1. Quick Test with Demonstration Data
By default, the master script is configured to run on the small-scale demonstration dataset provided in the `data/demonstration/` directory. This allows you to quickly verify that the pipeline executes correctly in your environment.

```bash
python Master_Pipeline_Integrated_Full_Strict_MultiBeta.py
```
> **Note:** The results generated from this demonstration data are random and do not reflect the actual biological findings reported in the paper. For a quick test, please ensure the SA parameters and grid search ranges are reduced as commented in the master script.

### 2. Reproducing the Study Results
To reproduce the exact results and figures reported in our study:

1. Download the full metagenomic datasets from the curatedMetagenomicData repository.
2. Place the formatted datasets into the `data/raw/` directory.
3. Modify the input data path in the master script to point to `data/raw/`.
4. Restore the original hyperparameter settings in the script (e.g., `SWEEPS_RANGE = ...`, `N_ENSEMBLE_TRIALS = 100`).
5. Run the script. The true results will be saved in the `Final_Results_SyNQA_Strict/` directory.

> **⚠️ Computational Cost Note:** > SyNQA relies on Simulated Annealing (SA) to solve an NP-hard combinatorial optimization problem (QUBO) and employs a strict Leave-One-Group-Out (LOGO) cross-validation from scratch to ensure unbiased feature selection. Therefore, reproducing the full study results with the original hyperparameter settings is highly computationally intensive and may take several hours (e.g., ~7 hours on an 11-core CPU).

## Output Details
When the pipeline completes, the following types of files are automatically generated in the `Final_Results_SyNQA_Strict/` directory:

* **Selected microbial taxa lists:** Details of the extracted microbial rewiring guild (e.g., `mechanism_detailed.csv`).
* **Performance metrics:** LOGO cross-validation evaluation results comparing SyNQA with baseline machine learning methods.
* **Network visualization figures:** Figures illustrating parameter landscapes, network properties, co-selection matrices, force-directed networks, and mechanistic correlations.

## Data Availability
The full metagenomic datasets analyzed in our study are publicly available in the [curatedMetagenomicData](https://waldronlab.io/curatedMetagenomicData/) repository (Bioconductor). To facilitate immediate testing and reproducibility, a script to generate a small-scale demonstration dataset is provided.

## Citation
If you use SyNQA in your research, please cite our paper:
> Arita, K., Nakano, Y., & Miyazaki, S. (2026). Synergistic Network QUBO Analysis (SyNQA): A combinatorial optimization framework for interaction-aware microbiome feature selection. *BMC Bioinformatics* (Under Review).

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
```
