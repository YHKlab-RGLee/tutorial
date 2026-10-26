# Local validation — 26 October 2026

The complete polymer workflow was executed successfully on Linux CPU using:
- Python 3.12.14, PyTorch 2.14.1+cpu
- NumPy 2.3.5, h5py 3.16.0, sisl 0.20.0
- GNU Fortran 13.3.0, serial SIESTA from commit fe3c3c18052b51101e438dad2ffa06f94e77b5bb
- System BLAS/LAPACK. The local compiler was unpacked in the workspace; explicit library paths and a compiler search path were used for linking. Colab uses the repository's standard build template.

The Colab runtime and its isolated Python 3.7.15 / PyTorch 1.12.1 environment have **not** been executed in this validation.

| Check | Result |
|---|---|
| Geometry | C8H16, 24 atoms; 20 x 20 x 10 Angstrom |
| Grid | 80 x 80 x 40 |
| SCF iterations | 9 |
| non-SCF iterations | 1 |
| Integrated SCF charge | 48.002403 electrons |
| Atomic density error | 17.301023% |
| ML density error | 0.820172% |
| SCF total energy | -1558.72530316 eV |
| non-SCF total energy | -1558.72283752 eV |
| Absolute energy error per atom | 0.102735 meV/atom |
| SCF sampled Gamma-X gap | 7.853300 eV |
| non-SCF sampled Gamma-X gap | 7.854800 eV |

Checked preprocessing output, matching feature scale, finite model inputs/outputs, RHO header and binary round trip, reference convergence, one-step completion, and the resulting density/band figures. Structure visualization was also constructed with sisl. Values may vary across environments. These results are for this example and do not establish basis, mesh, k-point convergence, force accuracy, or general model accuracy.
