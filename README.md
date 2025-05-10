## SINC: Standardization and Isotope Normalization for CSIA (with integrated correction and uncertainty quantification)

SINC is an open-source MATLAB package for processing compound-specific isotope analysis (CSIA) data from GC-IRMS. It provides a standardized workflow that includes drift correction, size normalization, scale correction, and uncertainty propagation, ensuring reproducibility and comparability across datasets.

This repository contains the MATLAB functions, tutorial, and example data used to apply SINC to _n_-alkane isotope measurements. The approach behind SINC is described in detail in:

📄 [Paper in review]

For a step-by-step tutorial, see:
📖 SINC_tutorial.m

## Key Features & Methodology
SINC implements several improvements over traditional CSIA-IRMS processing workflows. The major components include:

1️⃣ Drift Correction

Two options available:
- Compound-Specific Drift Correction ('CompoundSpecific') (Default)
- Global Drift Correction ('GlobalDrift')

Why It Matters: Traditional drift correction assumes all compounds drift uniformly, but in reality, different analytes can exhibit different drift patterns. SINC allows individualized drift corrections per compound, reducing systematic bias.

2️⃣ Scale Correction

Two options available:
- Regression-Based Scale Correction ('Regression') (Default)
- Two-Point Scale Correction ('TwoPoint')

Why It Matters: Traditional two-point scaling methods can introduce instability when reference standards exhibit high variability. SINC’s regression-based approach uses all available standard compounds, ensuring a more stable and accurate correction.

3️⃣ Outlier Detection & Handling

Configurable methods for outlier removal before correction:
- 'mad' (Median Absolute Deviation, default)
- 'zscore' (Removes values >3 standard deviations)
- 'iqr' (Interquartile Range)
- 'grubbs' (Detects single extreme outliers)
- 'none' (No outlier removal)

Why It Matters: Unfiltered extreme values can distort correction factors. SINC allows for flexible outlier detection methods to improve data quality.

4️⃣ Uncertainty Propagation

SINC propagates uncertainties from multiple sources:
Instrumental precision (repeatability)
- Drift correction uncertainty
- Size correction uncertainty
- Scale correction uncertainty
- Sample measurement uncertainty
- Reference gas uncertainty
  
Why It Matters: Many CSIA-IRMS studies underestimate uncertainty by considering only instrumental error. SINC integrates correction-based errors, producing more accurate confidence estimates.

5️⃣ Parameter Logging for Reproducibility

When enabled ('LogParams', 'yes'), SINC saves all input parameters to a timestamped log file (output/logs/).

Why It Matters: Ensures full reproducibility of processing settings, useful for method comparisons, debugging, and publication records.

## How to Get Started

1️⃣ Download SINC
Clone this repository or download the package manually.

2️⃣ Follow the tutorial
See 📖 SINC_tutorial.m for detailed instructions.

3️⃣ Run the workflow
Use infoQ to set parameters and processQ to run the pipeline.

🔗 Citation & References
If you use SINC in your research, please cite:
📄 [Paper in review]

For more details on the methodology and validation, refer to the manuscript.

👥 Contributors & Contact
SINC is developed and maintained by Julian Traphagan. If you have any questions, issues, or suggestions, feel free to open an issue on GitHub or reach out at jtraph1@lsu.edu.
