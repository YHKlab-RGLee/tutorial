# DeepSCF polymer tutorial inputs

The first notebook cell clones this repository into `/content/tutorial`.
Inputs and helpers are read directly from `/content/tutorial/deepscf/01/src`.
Run the accompanying notebook in order. All explanations are in the notebook.

## Files
- `inputs/STRUCT.fdf`: unchanged repository `polymer.fdf` (C8H16, 24 atoms).
- `inputs/C.psf`, `inputs/H.psf`: unchanged repository pseudopotentials.
- `inputs/common.fdf`: shared DZP/PBE settings, 1x1x8 k mesh, 80x80x40 real-space mesh, and Gamma-X bands.
- `inputs/scf/RUN.fdf`: reference SCF calculation.
- `inputs/non_scf/RUN.fdf`: one-step calculation reading `In.RHO`.
- `run_siesta.py`: execution and completion checks.
- `predict_polymer.py`: CPU inference with repository ResUNet, (R+S)x3 weights/configuration, and the reference grid header.
- `analysis.py`: density projection, full-grid density error, band plots, and total-energy comparison.

The polymer is periodic along z. X is explicitly defined as (0,0,1/2) in reciprocal-lattice coordinates; the supplied 10 Angstrom repeat cell produces folded bands. The 0.25 Angstrom spacing follows the pretrained examples, not a mesh-convergence study. MeshSizes overrides the requested MeshCutoff. Basis and k-point convergence are also not established.

Only feature enters the network; target is used for evaluating errors. The current notebook evaluates gaps along the sampled Gamma-X path, not a separate full-BZ gap search. All plotted bands share the SCF VBM as energy zero; sisl's per-file Fermi shift is undone first. The non-SCF output density is not used as the ML prediction.

The model archive is downloaded at runtime. No weights, compiled executables, or calculated grids are included. For a new attempt, rerun the preparation cell to create a fresh working directory, then rerun subsequent steps.

## Sources
Repository snapshot: `fe3c3c18052b51101e438dad2ffa06f94e77b5bb`.
- https://github.com/yhkimlab/DeepSCF/blob/fe3c3c18052b51101e438dad2ffa06f94e77b5bb/tutorial/1.preprocess/structures/examples/polymer.fdf
- https://github.com/yhkimlab/DeepSCF/tree/fe3c3c18052b51101e438dad2ffa06f94e77b5bb/tutorial/1.preprocess/input
- https://github.com/yhkimlab/DeepSCF/blob/fe3c3c18052b51101e438dad2ffa06f94e77b5bb/configs/predict/examples/r_plus_s_x3.yaml
- https://sisl.readthedocs.io/en/latest/api/io/generated/sisl.io.siesta.bandsSileSiesta.html

Repository license is retained as `LICENSE-DeepSCF`.
