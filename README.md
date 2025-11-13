
---

# Beyond Farm Size: Spatial Determinants of Cocoa Productivity in Ghana's Ashanti Region

[![License](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![R](https://img.shields.io/badge/R-4.3.1+-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)

## Overview

This repository contains the complete analytical pipeline for a comprehensive spatial econometric study examining cocoa productivity patterns across 2,612 georeferenced farms in Ghana's Ashanti Region. The research integrates GIS-based spatial analysis with advanced econometric modeling to identify determinants of yield variation in one of the world's most important cocoa-producing regions.

The project includes both **R-based spatial econometric analysis** (academic research) and a **Python-based interactive dashboard system** (practical implementation).

## Key Findings

- **Strong Spatial Clustering**: Cocoa yields exhibit pronounced spatial autocorrelation (Moran's I = 0.594, p < 0.001), with seven distinct productivity zones identified
- **Inverse Farm Size-Productivity Relationship**: Smallholdings (<2 ha) achieve 40-60% higher yields per hectare than larger operations (>5 ha)
- **Super-Elastic Tree Density Response**: A 10% increase in tree density translates to 12.9% yield improvement (elasticity = 1.29)
- **Two-Stage Production Process**: Strong spatial spillovers govern both production participation (ρ = 0.698) and conditional productivity (λ = 0.690)

## Interactive Dashboard

Explore the spatial patterns of cocoa productivity through our interactive web applications:

🌍 **[View R-based Interactive Dashboard](https://bwelson.github.io/spatialproject/)**

🐍 **[Launch Python Analysis System](#python-analysis-system)** (Run locally with your data)

## Repository Structure

```
├── README.md                                    # This file
│
├── R_Analysis/                                  # Academic Research Pipeline
│   ├── spatial_analysis_script.R
│   ├── model_comparison_script.R
│   ├── influential_observations_analysis.R
│   ├── point_pattern_analysis.R
│   
│
├── Python_Analysis/                             # Interactive Dashboard System
│   ├── perfect_cocoa_analyzer.py                # Main analysis system
│  
│
├── docs/
│   └── Beyond_Farm_Size_Spatial_Determinants.pdf  # Published manuscript

## Dual Analysis Framework

### 1. R-Based Spatial Econometric Analysis (Academic)

**Purpose**: Rigorous statistical modeling for academic publication

**Key Features**:
- Two-hurdle spatial regression framework
- Spatial autoregressive probit (SAR Probit)
- Spatial error model (SEM)
- Sensitivity analysis
- Point pattern analysis (Ripley's K, G-function)

**Technologies**: R 4.3.1, spdep, spatialreg, spatialprobit, spatstat

### 2. Python-Based Interactive System (Practical)

**Purpose**: User-friendly analysis and visualization for practitioners

**Key Features**:
- Automatic file format detection (CSV, Excel)
- Intelligent column mapping
- Comprehensive data cleaning
- Interactive folium maps
- Machine learning predictions
- Actionable recommendations
- Executive summary generation

**Technologies**: Python 3.8+, GeoPandas, Scikit-learn, Folium, Plotly

---

## Installation

### R Environment Setup

#### Required R Packages

```r
# Spatial analysis
install.packages(c("sf", "spdep", "spatialreg", "spatialprobit", "spatstat"))

# Visualization
install.packages(c("ggplot2", "tmap", "ggmap"))

# Data manipulation
install.packages(c("dplyr", "tidyr", "broom"))

# Model diagnostics
install.packages(c("car", "lmtest", "texreg"))
```

#### System Requirements

- R version 4.3.1 or higher
- Minimum 8 GB RAM (16 GB recommended for MCMC estimation)
- GDAL/PROJ libraries for spatial data processing

### Python Environment Setup

#### Using pip

```bash
# Create virtual environment
python -m venv cocoa_env
source cocoa_env/bin/activate  # On Windows: cocoa_env\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

#### Using conda

```bash
# Create conda environment
conda create -n cocoa_analysis python=3.8
conda activate cocoa_analysis

# Install packages
conda install geopandas pandas numpy matplotlib seaborn scikit-learn
pip install folium plotly openpyxl
```

#### requirements.txt

```
pandas>=1.3.0
numpy>=1.21.0
geopandas>=0.10.0
matplotlib>=3.4.0
seaborn>=0.11.0
folium>=0.12.0
plotly>=5.0.0
scikit-learn>=1.0.0
scipy>=1.7.0
shapely>=1.8.0
networkx>=2.6.0
openpyxl>=3.0.0
xlrd>=2.0.0
```

---

## Usage

### R-Based Spatial Analysis

#### 1. Load and Prepare Data

```r
source("R_Analysis/spatial_analysis_script.R")

# Read georeferenced data
cocoa_sf <- st_read("data/cocoa.shp")

# Check CRS
st_crs(cocoa_sf)
```

#### 2. Spatial Autocorrelation Analysis

```r
# Global Moran's I
moran_test <- moran.test(cocoa_sf$yld_pr_h, spatial_weights)

# Local Moran's I (LISA)
lisa_results <- localmoran(cocoa_sf$yld_pr_h, spatial_weights)
```

#### 3. Point Pattern Analysis

```r
# Ripley's K-function
K_function <- envelope(point_pattern, fun = Kest, nsim = 199)

# Nearest neighbor G-function
G_function <- envelope(point_pattern, fun = Gest, nsim = 199)
```

#### 4. Two-Hurdle Spatial Model

```r
# Hurdle 1: Spatial Probit (Participation)
sar_probit <- sarprobit(
  formula = yield_binary ~ farm_size + tree_density + ...,
  data = cocoa_sf,
  W = spatial_weights_matrix,
  ndraw = 5000
)

# Hurdle 2: Spatial Error Model (Conditional Yield)
spatial_error <- errorsarlm(
  formula = log(yield) ~ farm_size + tree_density + ...,
  data = cocoa_positive,
  listw = spatial_weights,
  method = "eigen"
)
```

#### 5. Model Comparison

```r
source("R_Analysis/model_comparison_script.R")

# Compare OLS, SAR, and SER
model_comparison <- compare_spatial_models(
  ols_model, sar_model, ser_model,
  metrics = c("AIC", "BIC", "RMSE", "Moran_I")
)
```

### Python-Based Interactive Analysis

#### 1. Basic Usage

```python
from perfect_cocoa_analyzer import PerfectCocoaFarmAnalyzer

# Initialize analyzer (supports CSV and Excel)
analyzer = PerfectCocoaFarmAnalyzer("data/cocoa_farms.xlsx")

# Run complete analysis pipeline
results = analyzer.run_perfect_analysis()
```

#### 2. Step-by-Step Analysis

```python
# Load and clean data
analyzer.detect_file_type_and_load()
analyzer.map_columns()
analyzer.comprehensive_data_cleaning()

# Create spatial features
analyzer.create_geodataframe()

# Perform analyses
analyzer.advanced_exploratory_analysis()
analyzer.intelligent_clustering_analysis(max_clusters=8)
analyzer.hotspot_analysis()

# Build predictive models
analyzer.predictive_modeling_suite()

# Generate visualizations
analyzer.create_stunning_visualizations()
analyzer.create_interactive_dashboard_map()

# Generate recommendations
analyzer.generate_actionable_recommendations()
analyzer.create_executive_summary_report()

# Export everything
export_path = analyzer.export_comprehensive_results()
```

#### 3. Access Analysis Results

```python
# View exploratory analysis
exploratory_results = analyzer.analysis_results['exploratory_analysis']

# View clustering results
clusters = analyzer.analysis_results['clustering_analysis']

# View model performance
models = analyzer.analysis_results['predictive_models']

# View recommendations
recommendations = analyzer.analysis_results['recommendations']
```

#### 4. Jupyter Notebook Example

```python
import pandas as pd
from perfect_cocoa_analyzer import PerfectCocoaFarmAnalyzer

# Load data
analyzer = PerfectCocoaFarmAnalyzer("data/cocoa_data.csv")

# Run analysis
results = analyzer.run_perfect_analysis()

# Display interactive map in notebook
analyzer.analysis_results['interactive_map']

# View recommendations
pd.DataFrame(
    analyzer.analysis_results['recommendations']['immediate_actions']
)
```

---

## Python Analysis System Features

### 🎯 Smart Data Processing

- **Automatic Format Detection**: Handles CSV, Excel (.xlsx, .xls)
- **Intelligent Column Mapping**: Automatically detects GPS coordinates, yield, farm size, etc.
- **Robust Data Cleaning**: 
  - GPS validation for Ghana boundaries
  - Outlier detection and removal
  - Missing value imputation
  - Data quality reporting

### 📊 Comprehensive Analysis

1. **Exploratory Analysis**
   - District-level statistics
   - Gender gap analysis
   - Farm size distribution
   - Correlation matrices

2. **Spatial Clustering**
   - K-means with optimal cluster selection
   - Silhouette score optimization
   - Cluster characterization

3. **Hotspot Identification**
   - Productivity quartile analysis
   - Geographic hotspot mapping
   - Performance comparison

4. **Predictive Modeling**
   - Random Forest regression (volume, yield)
   - Classification (productivity levels)
   - Feature importance analysis

### 🗺️ Interactive Visualizations

- **Multi-layer Folium Maps**:
  - Productivity levels (color-coded)
  - Spatial clusters
  - Hotspot markers
  - Production heatmaps
  - Multiple basemaps (OpenStreetMap, Terrain, etc.)

- **Static Visualizations**:
  - 16-subplot comprehensive dashboard
  - District comparisons
  - Gender analysis
  - Correlation heatmaps
  - Model performance charts

### 💡 Actionable Recommendations

Automatically generates categorized recommendations:
- **Immediate Actions** (3-6 months)
- **Short-term Strategies** (6-12 months)
- **Long-term Initiatives** (2-5 years)
- **Policy Recommendations**

Each recommendation includes:
- Priority level
- Issue description
- Specific actions
- Resource requirements
- Timeline
- Expected impact
- Cost estimate
- KPIs for monitoring

### 📤 Comprehensive Exports

- Cleaned datasets (CSV, GeoJSON)
- Interactive HTML maps
- Analysis reports (JSON)
- Executive summaries
- Model performance metrics
- District/gender analysis tables
- Recommendations documents

---

## Data Description

### Sample Characteristics

| Variable | Mean | SD | Min | Max |
|----------|------|-----|-----|-----|
| Farm size (ha) | 8.45 | 5.13 | 0.33 | 27.13 |
| Tree density (trees/ha) | 5.26 | 11.91 | 0 | 358.21 |
| Yield per hectare (kg/ha) | 40.45 | 63.77 | 0 | 233.92 |
| Household size (persons) | 5.20 | 2.07 | 1 | 13 |

### Zero-Inflation Structure

- **Non-producing farms** (yield = 0): 1,434 (54.9%)
- **Producing farms** (yield > 0): 1,178 (45.1%)
  - Mean yield: 89.68 kg/ha
  - Median yield: 73.72 kg/ha

### Required Data Format

#### For R Analysis
Shapefile (.shp) or CSV with columns:
- `latitude`, `longitude` (decimal degrees)
- `yld_pr_h` (yield per hectare)
- `farm_sz` (farm size in hectares)
- `tr_dnst` (tree density)
- Additional covariates as needed

#### For Python Analysis
CSV or Excel file with flexible column names:
- GPS coordinates (latitude/longitude, lat/lon, etc.)
- Production data (volume, yield, harvest_kg, etc.)
- Farm characteristics (size, trees, household size, etc.)
- Optional: district, gender, age, etc.

---

## Key Results

### 1. Spatial Dependencies

Both production stages exhibit strong spatial autocorrelation:

| Stage | Parameter | Estimate | p-value | Interpretation |
|-------|-----------|----------|---------|----------------|
| Participation | ρ (SAR) | 0.698 | < 0.001 | Strong neighbor spillovers |
| Productivity | λ (SEM) | 0.690 | < 0.001 | Spatially correlated unobservables |

### 2. Farm-Level Determinants (Hurdle 2)

| Variable | Coefficient | Elasticity | p-value |
|----------|-------------|------------|---------|
| log(Tree density) | 1.290 | 1.29 | < 0.001 |
| log(Yield per tree) | 0.975 | 0.97 | < 0.001 |
| Farm size | -0.029 | -2.9% per ha | < 0.001 |
| Previous FDP | -0.084 | -8.1% | 0.002 |

### 3. Model Performance

| Metric | Value |
|--------|-------|
| R² (original scale) | 0.927 |
| RMSE | 19.10 kg/ha |
| MAE | 7.08 kg/ha |
| Moran's I (residuals) | -0.016 (p = 0.858) |

### 4. Cluster-Level Productivity

| Cluster | N Farms | Mean Yield (kg/ha) | Dominant Districts |
|---------|---------|--------------------|--------------------|
| 5 | 259 | 171.29 | New Edubiase B/A |
| 2 | 295 | 127.85 | Asante Bekwai, Nsokote |
| 3 | 652 | 25.35 | Antoakrom A/B |
| 1 | 285 | 0.00 | Konongo A, Juaso |

---

## Policy Implications

### 1. Geographically Targeted Interventions

The identification of seven distinct productivity zones supports zone-specific strategies:

- **High-productivity clusters**: Focus on intensification (fertilizer, improved varieties)
- **Low-productivity clusters**: Prioritize foundational support (training, infrastructure)

### 2. Smallholder Prioritization

The inverse farm size-productivity relationship (β = -0.029, p < 0.001) provides empirical justification for targeting smallholders:

- Farms <2 ha achieve 40-60% higher yields per hectare
- Maximize productivity gains per unit investment

### 3. Replanting Programs

The super-elastic tree density effect (elasticity = 1.29) suggests substantial returns to intensification:

- 10% tree density increase → 12.9% yield improvement
- At current prices (GH¢56.64/kg), translates to GH¢655.33 additional revenue per hectare

### 4. Extension Service Redesign

Spatial spillovers (ρ = 0.698) indicate that interventions generate positive externalities:

- Target entire clusters rather than individual farms
- Leverage peer learning effects
- Design spatially coordinated programs

---

## Example Outputs

### Python System Generates:

```
perfect_cocoa_analysis/
├── export_20250113_143052/
│   ├── cleaned_cocoa_data.csv
│   ├── cocoa_farms_spatial.geojson
│   ├── interactive_dashboard_map.html
│   ├── recommendations.json
│   ├── executive_summary.json
│   ├── model_performance.json
│   ├── district_volume_stats.csv
│   ├── district_yield_stats.csv
│   ├── gender_analysis.json
│   ├── clustering_analysis.json
│   └── README.md
└── comprehensive_analysis.png
```

### Sample Recommendation Output:

```json
{
  "priority": "HIGH",
  "category": "District Support",
  "issue": "Low productivity in districts: Juaso, Konongo A",
  "recommendation": "Deploy agricultural extension officers to provide targeted training",
  "resources_needed": "Extension officers, training materials, soil testing kits",
  "timeline": "3-6 months",
  "expected_impact": "15-25% yield improvement",
  "cost_estimate": "Medium",
  "kpi": "Average yield per hectare in target districts"
}
```

---

## Reproducibility

### Computational Environment

**R Environment:**
```r
sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.3 LTS
```

**Python Environment:**
```python
import sys
print(f"Python {sys.version}")
# Python 3.8+ recommended
```

### Random Seeds

All stochastic procedures use fixed seeds for reproducibility:

**R:**
```r
set.seed(123)  # Spatial analysis
set.seed(456)  # MCMC estimation
set.seed(789)  # Bootstrap confidence intervals
```

**Python:**
```python
np.random.seed(42)  # NumPy operations
random_state=42     # Scikit-learn models
```

---

## Troubleshooting

### Common Issues

#### R Issues

**Problem**: Package installation fails
```r
# Solution: Install system dependencies (Ubuntu/Debian)
sudo apt-get install libgdal-dev libproj-dev libgeos-dev libudunits2-dev
```

**Problem**: MCMC takes too long
```r
# Solution: Reduce draws or use fewer neighbors
sarprobit(..., ndraw = 1000)  # Reduce from 5000
knn_obj <- knearneigh(coords, k = 5)  # Reduce from 7
```

#### Python Issues

**Problem**: GPS coordinates not detected
```python
# Solution: Check column names and manually map
analyzer.data.rename(columns={'lat': 'latitude', 'lng': 'longitude'}, inplace=True)
```

**Problem**: Memory error with large datasets
```python
# Solution: Sample data or increase memory
analyzer.data = analyzer.data.sample(n=1000)  # Reduce sample size
```

**Problem**: Missing dependencies
```bash
# Solution: Install specific package
pip install geopandas  # For spatial analysis
pip install openpyxl   # For Excel support
```

---

## Contributing

Contributions are welcome! Please follow these guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Commit your changes (`git commit -am 'Add new diagnostic'`)
4. Push to the branch (`git push origin feature/improvement`)
5. Open a Pull Request

### Development Setup

```bash
# Clone repository
git clone https://github.com/yourusername/cocoa-spatial-analysis.git
cd cocoa-spatial-analysis

# Install development dependencies
pip install -r requirements-dev.txt

# Run tests
pytest tests/

# Check code style
flake8 Python_Analysis/
black Python_Analysis/
```

---

## Citation

If you use this code or data in your research, please cite:

```bibtex
@article{bentum2025beyond,
  title={Beyond Farm Size: Spatial Determinants of Cocoa Productivity in Ghana's Ashanti Region},
  author={Bentum, Welson and Owiredu, Emmanuel Odame and Frimpong, Nana Kena and 
          Ofosu, Michael Asante and Boateng, Esther and Owusu, Kwaku Obeng and 
          Shadrack, Agyeman},
  journal={Journal Not Specified},
  volume={1},
  number={0},
  year={2025},
  publisher={Publisher Not Specified},
  doi={10.3390/1010000}
}
```

---

## License

This project is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0). See [LICENSE](https://creativecommons.org/licenses/by/4.0/) for details.

---

## Contact

**Welson Bentum**  
Department of Statistics and Actuarial Science  
Kwame Nkrumah University of Science and Technology, Kumasi, Ghana  
Email: bwelson523@gmail.com

---


## Related Resources

- [Spatial Econometrics Toolbox](https://www.spatial-econometrics.com/)
- [GeoPandas Documentation](https://geopandas.org/)
- [R Spatial Task View](https://cran.r-project.org/web/views/Spatial.html)
- [Ghana Cocoa Board](https://cocobod.gh/)

---

**Last Updated**: January 13, 2025  
**Status**: Published  
**Version**: 2.0.0 (Added Python analysis system)

---

⭐ **Star this repository** if you find it useful for your agricultural research or spatial analysis projects!
