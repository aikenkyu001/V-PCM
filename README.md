# V-PCM Engine: Virtual Photonic Computing Machine
**~ From Mathematical Axioms to Virtual Photonic Intelligence ~**

[![Academic Paper](https://img.shields.io/badge/Paper-V--PCM_en.md-blue.svg)](V-PCM_en.md)
[![Implementation](https://img.shields.io/badge/Implementation-Python%20%2F%20Fortran-green.svg)](Scripts/)
[![Theory](https://img.shields.io/badge/Theory-PKGF-orange.svg)](V-PCM_Plan.md)

## 🌌 Overview
**V-PCM (Virtual Photonic Computing Machine)** is a next-generation computing engine that realizes the **Parallel Key Geometric Flow (PKGF)** theory within a virtual software environment. Unlike traditional von Neumann architectures, V-PCM models intelligence as a continuous **Geometric Flow**, achieving $O(1)$ physical parallelism by emulating optical phenomena such as diffraction, interference, and dissipation.

**All experimental environments and source codes are publicly available at:**  
👉 [https://github.com/aikenkyu001/V-PCM](https://github.com/aikenkyu001/V-PCM)

This repository serves as a **Constructive Demonstration of Realizability** of PKGF theory, providing both a high-level Python prototype for visualization and a high-performance Fortran implementation for rigorous numerical verification.

## 🧠 Core Theory: PKGF
The engine is governed by the PKGF axiomatic system (**A1–U6**), which unifies three distinct phases of intelligence:
- **Constructive (C1–C3)**: Structure generation via semantic potentials.
- **Destructive (D1–D5)**: Information consolidation and entropy increase.
- **Unified (U1–U6)**: Spatiotemporal metabolism, spontaneous symmetry breaking, and dimensional jumps.

## 🚀 Key Features
- **Spatiotemporal Geometry**: Full implementation of spatial connection ($\nabla$), temporal connection ($\nabla_t$), and curvature ($F$).
- **O(1) Parallelism**: Emulation of photonic convolution using Gaussian kernels ($A$).
- **Dual-Language Verification**: Identical mathematical logic implemented in both Python (NumPy) and Fortran (LAPACK/Accelerate).
- **Optical Inverse Mapping**: Automated derivation of physical hardware parameters (Lens F-number, PSF $\sigma$) from virtual models.
- **Topological Analysis**: Tracking of structural components and spectral gaps during phase transitions.

## 📁 Repository Structure
```text
.
├── v_pcm.py                # Python Prototype (Visualization & ONP Analysis)
├── vpcm_fortran.f90        # High-Performance Fortran Implementation
├── V-PCM_jp.md             # Academic Paper (Japanese)
├── V-PCM_en.md             # Academic Paper (English - Native Level)
├── V-PCM_Plan.md           # Original Research & Implementation Plan
├── evolution_metrics.csv   # Numerical results (Python)
├── evolution_metrics_fortran.csv # Numerical results (Fortran)
├── summary.txt             # Experiment summaries and lens design
└──  exp*.png                # Visualization of metrics and matrix evolution
```

## 🛠 Installation & Usage

### Python Environment
Requirements: `python 3.x`, `numpy`, `matplotlib`
```bash
# Set up virtual environment
python3 -m venv venv
source venv/bin/activate
pip install numpy matplotlib

# Run simulation
python3 v_pcm.py
```

### Fortran Environment
Requirements: `gfortran`, `Accelerate Framework` (macOS) or `LAPACK`
```bash
# Build and Run
gfortran -O3 vpcm_fortran.f90 -framework Accelerate -o vpcm_fortran
./vpcm_fortran
```

## 📊 Experimental Results
The system demonstrates a clear **Geometric Phase Transition** at $t \approx 100$, where:
1. **Rank Jump (U6)** occurs as information organizes.
2. **Curvature ($F$) Diverges** during spontaneous symmetry breaking.
3. **Entropy ($S$) Saturates** at the theoretical maximum ($\ln N$).

![Matrix Evolution](exp1_matrix_evolution.png)
*Figure: Spatiotemporal self-organization of the Parallel Key ($K$).*

## 📜 Citation
If you use this engine or theory in your research, please cite:
> Miyata, F. (2026). Constructive Demonstration of Realizability and Numerical Verification of V-PCM: A Virtual Photonic Computing Engine Based on Parallel Key Geometric Flow (PKGF) Theory.

---
**Author**: Fumio Miyata  
