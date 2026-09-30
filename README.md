# Age-dependent association between hypertension and stroke risk in atrial fibrillation

Analysis code for a nationwide cohort of patients with atrial fibrillation. The outcome is ischemic stroke. Poisson regression on a Lexis data structure is used to estimate whether the association between hypertension and ischemic stroke changes with age. Age is modelled in bands and as a natural cubic spline. In one follow-up definition, oral anticoagulation is a time-varying state. In the other, follow-up ends at anticoagulation purchase.

## Files

- `HypertensionxAgeLL.R` — analysis script.

## Environment (R)

The cohort file is not part of this repository. Place `incident_cohort.sav` (one row per person, column names as in the script) in the working directory before running.

Install the CRAN packages loaded by the script (`Epi`, `haven`, `ggplot2`, `RColorBrewer`, `survival`, `knitr`, `gridExtra`). `splines` and `grid` ship with R. `cohorttools` is not on CRAN; install it before running. It supplies `mkratetable()`.

From this directory:

```bash
Rscript HypertensionxAgeLL.R
```

The script writes HTML rate and rate-ratio tables, PDF and JPEG figures, and `session_info.txt` into the working directory.

## License

MIT
