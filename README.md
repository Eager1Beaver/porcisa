# 🐖 Porcisa: Computational Porcine Heart–Torso Model for Early Ischemia

**Porcisa** is a modular C++/Python framework for simulating and visualizing electrophysiological dynamics of the **porcine heart-torso system** under ischemic conditions.  
It implements a cell-to-ECG multiscale model capturing the **triphasic conduction pattern** observed experimentally during early myocardial ischemia --- including **phase 1a conduction delay, interphase recovery,** and **phase 1b re-delay**.

---

## 🧩 Overview

This repository accompanies the paper:

> **Ilia Golub, Jan Azarov et al.**  
> *Computational Investigation of Early Ischemia: A Porcine Model Highlighting Activation Dynamics* (2025).  

The work introduces the first **species-specific porcine computational model** integrating:
- Ionic remodeling based on experimental porcine data  
- An anatomically realistic heart–torso mesh with rule-based fibre orientation  
- Bidomain simulations of electrical propagation using **Chaste**  
- Virtual 12-lead ECG computation  
- Python visualization of activation dynamics and conduction evolution  

---

## 🧱 Repository Structure

```
porcisa/
├── mesh/
│   ├── ecg_nodes.txt              # Virtual electrode node indices
│   └── pig_mesh_400k/             # Volumetric heart–torso mesh (not in repo)
│
├── output/                        # Simulation outputs (not tracked in Git)
│
├── paper/
│   ├── paper.pdf                  # Final manuscript
│   ├── figs/                      # Figures for the paper
│   └── latex_src/                 # LaTeX source files
│
├── scripts/                       # Python plotting utilities
│   ├── plot_convergence.py
│   ├── plot_singlecell.py
│   ├── prmtrs_at_evo.py
│   └── prmtrs_dr_evo.py
│
├── src/                           # C++ source modules
│   ├── PorcineConductivityModifier.hpp
│   ├── PorcineHeterogen.hpp
│   └── cell_factory/
│       └── PorcineCellFactory.hpp
│       └── ode/
│           ├── pig_ventr_ap_endo.cellml
│           ├── pig_ventr_ap_epi.cellml
│           └── pig_ventr_ap_m.cellml
│
├── test/                          # Unit and integration tests
│   ├── CMakeLists.txt
│   ├── ContinuousTestPack.txt     # Test list for batch runs
│   ├── TestSingleCell.hpp         # Steady-state AP & biomarker analysis
│   ├── TestConvergence.hpp        # Mesh resolution convergence
│   └── TestPorcineECG.hpp         # Whole-heart ischemic ECG simulation
│
├── .gitignore
├── porcisa.code-workspace
└── README.md
```

---

## ⚙️ Prerequisites

**Dependencies:**
- [Chaste](https://chaste.github.io/) (C++17, MPI-enabled build)
- CMake ≥ 3.12
- Boost ≥ 1.70
- PETSc + HDF5 support (for bidomain)
- Python ≥ 3.9 (for analysis scripts)
  - `numpy`, `matplotlib`, `pandas`

**Recommended:**
- Ubuntu 20.04 / 22.04  
- GCC ≥ 9.0  
- 32 GB RAM

---

## 🚀 Building and Running Tests

All tests are built within Chaste's CMake testing environment.

### 1. Configure
```bash
cmake ..
```

### 2. Compile a test
Example (single test):
```bash
make TestPorcineECG -j4
```

### 3. Run via `ctest`
```bash
ctest -V -R TestPorcineECG
```

---

### 🔄 Running Tests in Parallel

#### a. Parallel execution through `ctest`
```bash
cmake -D Chaste_NUM_CPUS_TEST=4 ..
make TestPorcineECG -j4
ctest -V -R TestPorcineECG
```

#### b. MPI execution
```bash
# Run on 1 core
./heart/test/TestPorcineECG

# Run on 4 cores
mpirun -np 4 ./heart/test/TestPorcineECG
```

---

## 🧪 Included Simulations

| Test | Description | Output |
|------|--------------|---------|
| **TestSingleCell.hpp** | Runs a single porcine ventricular cell to steady state; computes biomarkers such as APD, RMP, and peak voltage. | `output/singlecell/` |
| **TestConvergence.hpp** | Verifies activation time stability across different mesh resolutions in ischemic tissue. | `output/convergence/` |
| **TestPorcineECG.hpp** | Simulates activation and 12-lead ECGs across baseline, phase 1a, interphase, and phase 1b ischemic states. | `output/ecg/` |

---

## 📊 Visualization Scripts

Located in `/scripts/` (Python 3.x):

| Script | Purpose |
|--------|----------|
| `plot_singlecell.py` | Visualizes action potential traces and tabulates biomarkers. |
| `plot_convergence.py` | Plots activation-time vs. mesh-resolution curves. |
| `prmtrs_at_evo.py` | Illustrates evolution of activation time with ischemic parameters. |
| `prmtrs_dr_evo.py` | Plots dispersion of repolarization during ischemia. |

Run with:
```bash
python3 scripts/plot_singlecell.py
```

---

## 🧠 Scientific Summary

The **Porcisa** framework reproduces the experimentally observed *triphasic conduction pattern* during early ischemia:

1. **Phase 1a (1–10 min):** marked conduction delay and QRS widening driven by Na⁺ current suppression.  
2. **Interphase (10–20 min):** partial recovery due to transient balance among INa, INaL, and INa/K.  
3. **Phase 1b (20–30 min):** renewed slowing from rising tissue resistivity and current inhibition.

The model provides a mechanistic bridge between **ionic remodeling**, **tissue conductivity**, and **body-surface ECG morphology** in the porcine heart.

For detailed mathematical formulations and parameters, see **[`paper/paper.pdf`](paper/paper.pdf)** and **`Supplementary Materials.pdf`**.

---

<!--## 📚 Citation

If you use this repository, please cite:

> Golub I.A., Azarov J.E., Tsvetkova A.S., et al.  
> *Computational Investigation of Early Ischemia: A Porcine Model Highlighting Activation Dynamics.* (2025)  

---
-->
## 🧾 License

<!--© 2025 Ilia A. Golub and collaborators.-->  
This project is distributed for academic and non-commercial use under a **Creative Commons BY-NC-SA 4.0** license.

---

## 🐷 Acknowledgments

This work was supported by the **Program for Fundamental Research of the Russian Academy of Sciences (2022–2026)**, project # FUUU-2022-0068-1021052404529-3.

---
