# NekLBM

> **A High-Order Spectral Element Lattice Boltzmann Solver**

**NekLBM** is a high-order **Spectral Element Lattice Boltzmann Method (SELBM)** solver built upon the scalable **Nek5000** framework. By leveraging Nek5000’s robust mesh generation, high-performance I/O, and advanced MPI distribution systems, NekLBM provides a powerful platform for simulating complex fluid dynamics.

---

## 🚀 Key Capabilities

NekLBM is engineered for high-fidelity simulations across multiple regimes:

* **Single-phase flows**
* **Multi-phase flows**
* **Multi-species fluid dynamics**

## 🏗️ Framework Integration

NekLBM is fully integrated with the **Nek5000** ecosystem, inheriting its industrial-grade capabilities for High-Performance Computing (HPC):

* **Mesh Geometry:** Utilizes native spectral element mesh information for complex geometries.
* **I/O System:** High-efficiency data handling for large-scale simulations.
* **Distributed Systems:** Optimized for massive parallelism via MPI.



## 📚 Citation

If you use **NekLBM** in your research, please cite our latest work:

> [1] **Zhao, C., Patel, S.S., Lin, H.L., Min, M. and Lee, T., 2026.** A flux bounce-back scheme for the filtered Spectral Element Lattice Boltzmann Method. *Computers & Fluids*, p.106987. [https://doi.org/10.1016/j.compfluid.2025.106987](https://doi.org/10.1016/j.compfluid.2025.106987)

---

# 🛠️ Installation

**NekLBM** is built on top of the Nek5000 framework. Before proceeding, ensure you have a working environment for the core solver.

## 1. Prerequisites (Nek5000)
Follow the official **Nek5000 Quickstart Guide** to install the base framework, compilers, and dependencies:
👉 [Nek5000 Installation & Quickstart](https://nek5000.github.io/NekDoc/quickstart.html)

## 2. Core Modification
To enable the LBM driver, you must replace the default `drive1.f` file within the Nek5000 source tree:
* **Target Path:** `Nek5000/core/drive1.f`
* **Action:** Replace this file with the `drive1.f` provided in the NekLBM repository.

## 3. Lattice & Physics Configuration
The lattice structures and simulation parameters are defined in specific files provided in NekLBM/parameters. You have to copy those files to `Nek5000/core` folder. Choose the file corresponding to your simulation requirements:

### Single-Phase Lattices
| Dimension | Lattice Type | Parameter File |
| :--- | :--- | :--- |
| **2D** | D2Q9 | `LBMD2Q9` |
| **3D** | D3Q13 | `LBMD3Q13` |
| **3D** | D3Q15 | `LBMD3Q15` |
| **3D** | D3Q19 | `LBMD3Q19` |
| **3D** | D3Q27 | `LBMD3Q27` |

### Multi-Phase & Multi-Species
For specialized physics, parameters are stored in the following files:

* **Multi-Phase Flow:**
    * 2D: `LBMPHASE2D`
    * 3D: `LBMPHASE3D`
* **Multi-Species Flow:**
    * 2D: `LBMDIFFUSION2D`
    * 3D: `LBMDIFFUSION3D`
Once the files are exchanged and parameters are set:
1. Initialize your case directory using `makenek`.
2. Compile as you would for a standard Nek5000 project.
