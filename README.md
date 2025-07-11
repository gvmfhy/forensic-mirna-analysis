# Vibe Coding: AI-assisted Development of Forensic miRNA Analysis Pipeline using Claude Opus 4

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15860835.svg)](https://doi.org/10.5281/zenodo.15860835)

A bioinformatics pipeline for identifying microRNA (miRNA) signatures in forensic body fluid samples, developed through transparent human-AI collaboration. This project demonstrates "vibe coding" - using natural language to guide AI through complex scientific software development, resulting in a working pipeline that identifies 393 miRNA markers for forensic body fluid identification.

## 🎯 Development Transparency

This entire pipeline was developed through AI-assisted programming using Claude Opus 4. The complete, unedited development session is available:
- **[View the full vibe coding transcript](docs/development_logs/2025-07-11-session/vibe_coding_session.html)** - 338,952 lines of human-AI interaction
- **[Read about the development process](docs/development_logs/2025-07-11-session/README.md)** - Including failures and recovery strategies
- **[Visit the web version](https://rna-is-cool.netlify.app/)** - For rendered documentation

## 🔬 Project Overview

This pipeline processes two complementary microarray datasets:
- **GSE153135**: Two-color Agilent arrays (GPR format) with pooled samples
- **GSE49630**: Affymetrix ST arrays (CEL format) with individual samples

The analysis identifies 393 forensic miRNA marker candidates with large effect sizes suitable for body fluid identification.

## 📊 Key Results

From analysis of n=5 samples per body fluid type:
- **Blood markers**: 155 candidates (e.g., hsa-miR-486-5p with 3243-fold enrichment*)
- **Semen markers**: 167 candidates (e.g., hsa-miR-891a with 175-fold enrichment*)
- **Saliva markers**: 49 candidates (e.g., hsa-miR-205 with 275-fold enrichment*)
- **Vaginal markers**: 22 candidates (e.g., hsa-miR-138-1-star with 7-fold enrichment*)

*Note: Large fold-changes may reflect variance in small sample sizes and require validation

## ⚠️ Important Limitations

- **Small sample size**: Only n=5 samples per body fluid type
- **No independent validation**: Results require confirmation in larger cohorts
- **Statistical constraints**: FDR < 0.05 not achievable with 13,564 tests and small n
- **Forensic applicability**: These are research findings, not validated forensic markers
- **Development history**: See development logs for complete transparency including errors

## 🚀 Quick Start

### Prerequisites

- Python 3.8+ with pandas, numpy, scipy, matplotlib, seaborn
- R 4.5+ with BiocManager, oligo, limma packages
- macOS/Linux environment (tested on macOS ARM64)

### Installation

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/forensic-mirna-analysis.git
cd forensic-mirna-analysis

# Set up Python environment
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Set up R environment
Rscript setup_r_env.sh
```

### Data Download

Due to size constraints, raw data files are not included. Download from GEO:

```bash
# Create data directories
mkdir -p data/raw/GSE153135_GPR data/raw/GSE49630_CEL

# Download GSE153135 (GPR files)
# Visit https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153135

# Download GSE49630 (CEL files)  
# Visit https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49630
```

### Running the Analysis

```bash
# 1. Process GPR files (two-color arrays)
python scripts/preprocessing/gpr_parser.py

# 2. Process CEL files (Affymetrix arrays)
Rscript scripts/preprocessing/cel_processor_minimal.R

# 3. Run forensic marker analysis
python scripts/analysis/cel_practical_forensic.py
```

## 📁 Project Structure

```
├── data/
│   ├── raw/                    # Raw data files (not in repo)
│   └── processed/              # Normalized expression matrices
├── scripts/
│   ├── preprocessing/          # Data parsing and normalization
│   └── analysis/              # Statistical analysis and visualization
├── results/
│   ├── continuous_expression_analysis/  # GPR analysis results
│   └── cel_practical_forensic/         # CEL forensic markers
├── docs/
│   ├── DETAILED_ANALYSIS_WALKTHROUGH.md
│   └── RESULTS_AND_VISUALIZATIONS_STRUCTURE.md
└── requirements.txt
```

## 🔍 Analysis Approach

### Multi-Tier Expression Framework
The pipeline implements a multi-tier detection system for continuous expression data:
- **High confidence**: log2 ratio > -2.0
- **Moderate confidence**: -6.0 to -2.0  
- **Low confidence**: -10.0 to -6.0
- **Undetected**: < -10.0

### Statistical Methods
- Non-parametric Wilcoxon tests (suitable for small samples)
- Cohen's d effect size (prioritized over p-values)
- FDR correction for multiple testing
- Forensic-specific thresholds (>4-fold change, >99% specificity)

## 📚 Documentation

Comprehensive documentation is available:
- [Detailed Analysis Walkthrough](DETAILED_ANALYSIS_WALKTHROUGH.md) - Step-by-step process
- [Results Structure](RESULTS_AND_VISUALIZATIONS_STRUCTURE.md) - Output file descriptions
- Fully documented script versions in `*_documented.py/R` files

## ⚠️ Important Notes

### Platform-Specific Requirements
- **CEL files MUST use oligo package** (not affy) for ST arrays
- **macOS ARM64**: Requires specific R configuration (see setup_r_env.sh)

### Statistical Limitations
- Small sample sizes (n=2-5 per fluid) limit statistical power
- Strict FDR < 0.05 not achievable with 13,564 tests
- Focus on effect size over p-values for forensic relevance

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add comprehensive documentation
4. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see LICENSE file for details.

## 🙏 Acknowledgments

- GEO datasets GSE153135 and GSE49630 authors
- Bioconductor community for R packages
- Forensic miRNA research community
- Claude Opus 4 (Anthropic) for AI programming assistance

## 💡 LLMs and Scientific Discovery

Historically, experimental science has operated in siloed compartments, often insulated from broader technological innovations, especially those involving computational advancements. Computationally adjacent "dry-lab" fields rapidly incorporate cutting-edge AI tools, while traditional "wet-lab" research frequently lags behind due to habits of secrecy, selective reporting, and institutional inertia.

Despite known errors and current limitations, the widespread use of LLMs ensures their lasting impact. They underpin future technologies and scientific methodologies, making transparent practices today not just beneficial but essential.

Scientific publishing currently suffers from an entrenched culture of presenting only successful results. By openly documenting failures, misunderstandings, and friction points through transparent, raw development logs, we directly challenge this problematic status quo. Publishing these records online, such as full terminal logs hosted on publicly accessible websites, achieves two goals:

First, it immediately improves reproducibility, honesty, and transparency in scientific discourse. Researchers can directly examine genuine human-AI interactions, learn from recorded failures, and clearly understand breakdown points. We can better criticize each other and offer our domain specific knowledge when we notice errors.

Second, openly shared development logs represent an investment in future AI development. Public repositories containing candid transcripts may enter LLM training datasets, helping AI models learn from realistic, error-inclusive examples of human-AI collaboration, thereby improving their capability to navigate and mitigate such friction points.

By contributing realistic scientific problem-solving data to public repositories, individual researchers actively shape how future models understand and respond to scientific queries. This approach values instructive failures as essential components of scientific progress.

Publishing raw terminal logs directly targets the largest visibility gap in science by surfacing misfires and dead ends that traditional papers omit. This transparency provides peers—and future LLMs—access to the full causal chain of discovery, improving present-day scientific transparency and influencing human-AI interactions.

## 📞 Contact

For questions or collaborations, please open an issue on GitHub.

---

**Note**: This pipeline was developed for research purposes. Forensic applications require additional validation on independent samples. The development process demonstrates that sophisticated bioinformatics tools can be created through AI collaboration by researchers without traditional programming expertise.