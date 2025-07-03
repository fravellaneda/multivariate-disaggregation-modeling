# Multivariate Disaggregation Modeling

**Author**: Fernando Rodríguez Avellaneda, Erick A. Chacón-Montalván and Paula Moraga
**Paper**: Multivariate disaggregation modeling of air pollutants: a case‑study of PM2.5, PM10 and ozone prediction in Portugal and Italy ([arXiv:2503.12394](https://arxiv.org/abs/2503.12394))  

This repository implements the **simulation study**, accompanying the above paper. The goal is to demonstrate and validate the proposed **spatial multivariate disaggregation** approach using synthetic and real-world data.

---

## 📂 Repository Contents

```
- `create_dataset.R` — Simulates multivariate areal-level datasets
- `disaggregation_model.R` — Defines and fits the disaggregation model using INLA
- `aux_fun.R` — Auxiliary functions used to create or model the dataset
- `make_plots.R` — Generates visualizations for the simulation results
- `run_file.R` — Script to run the complete simulation pipeline
```

---

## 🚀 How to Run the Simulation Study

From the R console or terminal:

1. **Simulate synthetic data:**

   ```r
   source("create_dataset.R")
   ```

2. **Fit the model to the simulated data:**

   ```r
   source("disaggregation_model.R")
   ```

3. **Generate plots of the model output:**

   ```r
   source("make_plots.R")
   ```

Or run the entire pipeline via:

```r
source("run_file.R")
```

---

## 📝 Citation

If you use this code, please cite:

Rodríguez Avellaneda, F., Chacón‑Montalván, E. A., & Moraga, P. B. (2025).  
*Multivariate disaggregation modeling of air pollutants: a case‑study of PM2.5, PM10 and ozone prediction in Portugal and Italy*.  
arXiv:2503.12394

---

## 📬 Contact

**Fernando Rodríguez Avellaneda**  
📧 [fravellaneda@kaust.edu.sa](mailto:fernando.rodriguezavellaneda@kaust.edu.sa)

---
