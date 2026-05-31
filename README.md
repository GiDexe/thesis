# 🌍 ISOs Diffusion in the EU

**An Endogenous Network Recovery for Environmental and Quality Management Systems**

## 📖 Overview

Do voluntary standards diffuse across countries through a network of mutual
influence, and is that network the same for quality and for environmental
standards? This thesis recovers the unobserved diffusion network directly from
ISO adoption data across the 27 EU member states, for **ISO 9001** (Quality) and
**ISO 14001** (Environment), using the network recovery procedure of De Paula,
Rasul & Souza (2025) via the `recoverNetwork` package.

## 📂 Repository layout

```
MAIN.qmd                       Analysis + paper body
estimation.r                   Network recovery pipeline (cleaning + grid search)
paper.tex                      LaTeX paper, assembles tables/figures from assets/
thesis_helpers_refactored.r    Helper functions sourced by MAIN.qmd
Makefile                       Build targets
assets/                        Generated tables and figures
startup_scripts/               Dependency installer run on container start
.devcontainer/                 Dev Container definition
data/                          NOT in git — provided via Dropbox
```

## 💾 Data

The `data/` folder is **not** tracked by git. It is bind-mounted from a Dropbox
folder and contains the raw inputs, the cleaned panels, and the precomputed
estimates (`data/output/`). Since recovery runs ~72h, those estimates are the
canonical artefacts the paper is built from. Request the Dropbox link to obtain it.

## 🔁 Replication

### 1. ⚙️ Requirements

Developed on **Linux (Ubuntu)**. On **Windows**, use **WSL2** with Ubuntu and run
everything from the WSL shell. macOS works too (ARM64 supported). You need:

- 🐳 **Docker** (Engine on Linux/WSL2, or Docker Desktop)
- 🧩 **VS Code** + the **Dev Containers** extension
- 🌿 **git**
- 💾 the [**`data/`** folder](https://www.dropbox.com/scl/fo/wmldg06aek2m7xekbjr09/AG8_2z6wtD2aoNHxZwyc5BI?rlkey=ap8zxcg5ij4lp0mkgm5lzenlm&st=cmaov3l8&dl=0) (from Dropbox)
Everything else lives inside the container, built on
[**Academic Docker**](https://github.com/rferrali/AcademicDocker) 🎓 — a
reproducible image bundling R (rig + renv), Quarto, and TeXLive 2025.

### 2. 📥 Get the code and data

```bash
git clone https://github.com/GiDexe/thesis.git
cd thesis
# download the Dropbox `data` folder, then point it at ./data
```

### 3. 🔧 Configure the Dev Container

Edit `.devcontainer/devcontainer.json` and set the bind-mount paths (under WSL2,
use Linux paths like `/home/<user>/...`):

```jsonc
"mounts": [
  "source=devcontainer-renv-cache,target=/renv/cache,type=volume",
  "source=<ABSOLUTE PATH TO data>,target=/workspaces/thesis/data,type=bind",
  "source=<ABSOLUTE PATH TO startup_scripts>,target=/startup_scripts,type=bind"
]
```

### 4. 🚀 Open in the container

In VS Code: **Dev Containers: Reopen in Container**. On first launch the image is
pulled and the dependencies are installed automatically.

### 5. 🛠️ Build

```bash
make paper     # 📄 render MAIN.qmd, then compile paper.tex → paper.pdf
make results   # 📊 render MAIN.qmd only (assets + results.pdf)
make rerun     # 🔁 recompute estimates from raw data, then rebuild (~72h)
make clean     # 🧹 remove generated artefacts (never touches data/)
```

- **`make paper`** — normal entry point: builds the full PDF from the precomputed estimates.
- **`make results`** — refreshes tables/figures without recompiling the LaTeX paper.
- **`make rerun`** — recomputes the networks from scratch (needs raw `data/`; heavy).
- **`make clean`** — clears generated outputs only; leaves `data/` untouched.

## 🙏 Acknowledgements

Built on [rferrali/AcademicDocker](https://github.com/rferrali/AcademicDocker) 🎓.
Network recovery via the `recoverNetwork` implementation of De Paula, Rasul &
Souza (2025).
