# 🔭 Systematic Uncertainties in DES5yr & Pantheon+ SNe Ia and Their Impact on Cosmological Constraints

> **B.Sc. Thesis in Astronomy** — Universidad de Antioquia, 2026  
> **Author:** Melkyn Snneyder Quintana Rojas  
> **Supervisor:** Juan Carlos Muñoz-Cuartas, Ph.D.

---

## Overview

The **Hubble tension** — a ~5σ discrepancy between direct and indirect measurements of H₀ — is one of the biggest open problems in modern cosmology. This thesis investigates it by comparing two major Type Ia Supernovae datasets under a flat ΛCDM model:

| Dataset | SNe Ia | Redshift range | Light-curve model |
|---|---|---|---|
| **Pantheon+** (Riess et al. 2022) | 1,701 | 0.001 – 2.26 | SALT2 |
| **DES 5yr** (Abbott et al. 2024) | 1,829 | 0.01 – 1.12 | SALT3 |

The core question: **does the SALT2 → SALT3 transition introduce systematic differences in cosmological parameter estimation?**

---

## Methods

- MCMC sampling via [`emcee`](https://emcee.readthedocs.io/) with full covariance matrix (stat + sys uncertainties)
- Classical H₀ estimation with redshift binning + bootstrapping (Anderson-Darling test)
- SALT parameter comparison (c, x₀, x₁, mB) for common supernovae across datasets
- Host-galaxy cluster systematic: SDSS DR7 cross-match via `BallTree`

**Stack:** `Python · emcee · pandas · NumPy · SciPy · Matplotlib · Astropy`

---

## Key Results

| Dataset | H₀ (km/s/Mpc) | Ωₘ |
|---|---|---|
| Pantheon+ | 72.63 ± 0.56 | 0.3668 ± 0.0416 |
| Pantheon+ & Cepheids | 72.63 ± 0.56 | 0.3681 ± 0.0416 |
| DES 5yr | 69.54 ± 0.54 | 0.2762 ± 0.0222 |

The **color parameter (c)** shows a clear systematic offset between datasets for common supernovae, large enough to shift the distance modulus µ and propagate into cosmological inference. The H₀ difference persists across all analysis approaches (p < 0.05, Anderson-Darling).

---

## Repository Structure

```
├── data/
│   ├── Pantheon+/          # Pantheon+ dataset
│   ├── DES/                # DES 5yr dataset
│   └── Clusters/                # DES 5yr dataset
├── notebooks/
│   ├── 01_data_overview.ipynb
│   ├── 02_survey_analysis.ipynb
│   ├── 03_zbinning.ipynb
│   ├── 04_host_galaxy_clusters.ipynb
│   └── 05_salt_comparison.ipynb
├── src/
│   ├── cosmo_functions.py           # ΛCDM model & luminosity distance
│   ├── mcmc_functions.py          # Likelihood, Covariance regularization and emcee setup
│   ├── data_import.py                # Data reader
│   ├── crossmatching.py                # Crossmatch SN - Clusters
│   └── plot_functions.py               # Plotting
├── figures/                   # All plots and corner plots
├── Thesis_MelkynQuintana.pdf
└── README.md
```

---

## Key References

- Riess et al. (2022) — Pantheon+ dataset & H₀ = 73.04 ± 1.04 km/s/Mpc
- Abbott et al. (2024) — DES 5yr dataset
- Kenworthy et al. (2021) — SALT3 light-curve model
- Foreman-Mackey et al. (2013) — `emcee` MCMC library
- Muñoz-Cuartas et al. (2011) — SDSS DR7 galaxy cluster catalog
- Aghanim et al. (2020) — Planck 2018 cosmological parameters

---

## License

[Creative Commons BY-NC-SA](https://creativecommons.org/licenses/by-nc-sa/4.0/) — Universidad de Antioquia.