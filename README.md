<div align="center">

# 🌡️ Expired by Climate

### *A 45-Year Global Audit of Pharmaceutical Stability Test Conditions Against a Warming Climate*

[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Data: ERA5 + NOAA](https://img.shields.io/badge/Data-ERA5%20%2B%20NOAA-blue)](https://cds.climate.copernicus.eu/)
[![EurekaDev 2026](https://img.shields.io/badge/EurekaDev-2026%20Research%20Track-orange)](https://eurekadev2026.devpost.com)

---

**Every medicine you take has a shelf life — and that shelf life was proven inside a test chamber set to a temperature chosen decades ago. We checked whether those temperatures still hold.**

</div>

---

## 📌 The Problem

Before a medicine is approved, batches are stored for a year or more at a fixed "long-term" condition to prove the shelf life printed on the pack. The WHO sets three such conditions:

| Climatic Zone | Test Condition | Intended Coverage |
|:---:|:---:|:---|
| **Zone II** | 25 °C / 60 % RH | Temperate markets (US, EU, Japan) |
| **Zone IVa** | 30 °C / 65 % RH | Hot & humid regions |
| **Zone IVb** | 30 °C / 75 % RH | Hot & very humid regions |

These numbers were derived from **1980s–2000s climate data** using the Mean Kinetic Temperature (MKT) formula (Haynes 1971), with a safety margin added (Grimm 1998).

**The climate has warmed since.** If a region's effective temperature now exceeds its test condition, every shelf life approved there is *optimistic* — and degraded medicines under-dose patients.

---

## 🔬 What We Did

We treated each regulatory test condition as a **falsifiable hypothesis about the climate** and tested it with the regulators' own method — on current data.

For **18 places worldwide** (10 from ERA5 reanalysis, 8 from NOAA weather stations), we:

1. Computed **MKT and water-vapour partial pressure** for every year from **1981 to 2025**
2. Compared two 20-year windows: **1981–2000** vs. **2006–2025**
3. Measured the **safety margin** remaining between the actual climate and the test condition
4. Translated exceedances into **effective shelf-life loss** via the Arrhenius equation

All hypotheses and rejection rules were **pre-registered before the data run**.

---

## 📊 Key Findings

```
┌──────────────────────────────────────────────────────────────────────┐
│  ✦  Effective temperature ROSE in all 18 places (0.35–1.33 °C)     │
│  ✦  Every 95% bootstrap CI excludes zero                           │
│  ✦  7 of 18 places now EXCEED their test temperature               │
│  ✦  Positive control (Tokyo) passed — matching Miura 2024          │
└──────────────────────────────────────────────────────────────────────┘
```

| City | Zone | Test (°C) | MKT 2006–2025 (°C) | Margin (°C) | 24-mo label → |
|:---|:---:|:---:|:---:|:---:|:---:|
| **Dubai** | IV | 30.0 | 31.81 | **−1.81** | ~20 months |
| **Phoenix** | II | 25.0 | 29.35 | **−4.35** | ~15 months |
| **Houston** | II | 25.0 | 25.41 | **−0.41** | ~23 months |
| **Tampa** | II | 25.0 | 25.67 | **−0.67** | ~22 months |
| **Miami** | II | 25.0 | 26.60 | **−1.60** | ~21 months |
| Tokyo ✓ | II | 25.0 | 22.09 | +2.91 | 24 months ✓ |

> *After bias correction of +0.3 °C. Shelf-life estimates use default activation energy (83.1 kJ/mol).*

**On humidity:** Chennai, Colombo, and Yangon already hold more water vapour than the 30 °C / 65 % test chamber produces.

---

## 🏗️ Repository Structure

```
Expired-by-Climate/
│
├── code/
│   ├── analyse.py                    # Reproduces every number & figure in the paper
│   └── stability_drift_pipeline.py   # Gridded ERA5 pipeline (download → analyse)
│
├── data/
│   └── expired_by_climate_results.csv  # One row per place per year (MKT, Pd, source)
│
├── figures/
│   ├── fig1_model.png                # Conceptual model diagram
│   ├── fig2_margins.png              # Temperature margin dumbbell chart
│   ├── fig3_humidity.png             # Humidity margin chart (Zone IV)
│   └── fig4_shelflife.png            # Effective shelf-life chart
│
├── video/
│   ├── Expired_by_Climate_Video.mp4  # 3:52 presentation (1080p, captions burned in)
│   ├── Expired_by_Climate_Slides.pptx
│   └── Expired_by_Climate_Captions.srt
│
├── Expired_by_Climate_Research_Paper.pdf   # Full research documentation
├── Project_Description.txt                 # Devpost project description
└── README.md                              # ← You are here
```

---

## ⚡ Quickstart

### Prerequisites
```bash
pip install numpy pandas matplotlib
```

### Reproduce Every Result
```bash
# Runs analysis on the pre-computed CSV, outputs all figures + results.json
python code/analyse.py
```

### Run the Full ERA5 Pipeline (optional — requires CDS account)
```bash
pip install cdsapi xarray netCDF4

# Verify the MKT implementation passes all self-tests
python code/stability_drift_pipeline.py selftest

# Download ERA5 data for a region (example: India bounding box)
python code/stability_drift_pipeline.py download 1980 2025 --area 38 66 6 98

# Analyse downloaded NetCDF files
python code/stability_drift_pipeline.py analyse era5_hourofday_*.nc --out results.nc
```

---

## 🔧 Methods at a Glance

| Component | Method | Reference |
|:---|:---|:---|
| **Effective temperature** | Mean Kinetic Temperature (MKT) | Haynes 1971 |
| **Indoor floor** | Grimm 19 °C floor (heated storage) | Grimm 1998 |
| **Humidity** | Buck (1981) saturation vapour pressure | Buck 1981 |
| **Climate data** | ERA5 reanalysis + NOAA GHCN-Daily | Hersbach+ 2020, Menne+ 2012 |
| **Trend estimation** | Theil–Sen slope (robust to outliers) | — |
| **Uncertainty** | 5,000-iteration bootstrap CIs | — |
| **Shelf-life impact** | Arrhenius zero-order approximation | USP / ICH convention |
| **Bias correction** | +0.3 °C (daily max/min vs hourly) | Measured in this study |

---

## ✅ Built-in Validation

| Check | Purpose | Result |
|:---|:---|:---|
| **Positive control** | Tokyo should pass (Miura 2024 found no issue) | ✅ Passed (22.09 °C vs 25 °C) |
| **Bias measurement** | Daily max/min MKT vs hourly station MKT | +0.3 °C bias → corrected |
| **ERA5 vs Station** | Cross-validate reanalysis against ground truth | Agreement within expected range |
| **ICH worked example** | MKT of [20, 40] °C should ≈ 34.4 °C | ✅ 34.4 °C (selftest) |

---

## 📝 Citation

If you use this data or code in your work, please cite:

```
Expired by Climate: A 45-Year Global Audit of Pharmaceutical Stability
Test Conditions Against a Warming Climate. EurekaDev 2026 Research Track.
https://github.com/AdityaChandel11/Devpost-Eureka2026
```

---

## 📚 Key References

- Grimm W. (1998) *Drug Dev Ind Pharm* 24(4):313–325. [doi:10.3109/03639049809085626](https://doi.org/10.3109/03639049809085626)
- Haynes JD. (1971) *J Pharm Sci* 60:927–929. [doi:10.1002/jps.2600600629](https://doi.org/10.1002/jps.2600600629)
- Hersbach H, et al. (2020) *QJRMS* 146:1999–2049. [doi:10.1002/qj.3803](https://doi.org/10.1002/qj.3803)
- Kaplan WA, Hamer DH, Shioda K. (2025) *One Health* 20:100957. [doi:10.1016/j.onehlt.2024.100957](https://doi.org/10.1016/j.onehlt.2024.100957)
- Miura E. (2024) *Ther Innov Regul Sci* 58:184–191. [doi:10.1007/s43441-023-00584-4](https://doi.org/10.1007/s43441-023-00584-4)
- WHO (2018) TRS 1010, Annex 10: Stability testing of APIs and finished pharmaceutical products.

---

<div align="center">

*A standard is a forecast. Forecasts need checking.*

**Made for [EurekaDev 2026](https://eurekadev2026.devpost.com) · Research Track · Biology/Medical and Environmental Science**

</div>
