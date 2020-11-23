# Introduction

This is a collection of code that was used for the Ph.D. thesis

    Computational Gene Expression Deconvolution

by Dominik J. Otto. This repo only contains the code used in the example
application "inference of cell decomposition". The code that was applied
applied for the deoncolution of expression patterns is collected in another
repo.

# What to run

Make sure a directory `../data` is available.
Run everything in `main.sh` to make the model and `docker/Readme` to make
and upload the docker images.

# Dependencies

We us the [Slurm Workload Manager](https://slurm.schedmd.com/documentation.html)
to shechule jobs on our cluster system with `sbatch`. However, all scriptes
called with `sbatch <script name>` can be run without `sbatch`.

To make R 3.5.2 available we use the
[Lmod module system](https://lmod.readthedocs.io/en/latest/).
The respective calls `source /etc/profile.d/modules.sh`
and `module load R/3.5.2-0` can be omitted if R 3.5.2
is already available.

The calls `source ~dominik.otto/enterEnv.sh` load a
specially build
[conda environment](https://docs.conda.io/projects/conda/en/latest/).
If all use python packages are availabe these calls can be omitted.

The python multiprocessing package had a bug not allowing to share
files larger than 4GB between the threads. Please make sure
this bug is fiexed, e.g., by implementing this PR in your version
of multiprocessing: https://github.com/python/cpython/pull/10305

# License

Copyright (C) 2019 Fraunhofer-Gesellschaft zur Foerderung der angewandten
Forschung e.V. acting on behalf of its Fraunhofer Institute for Cell Therapy
and Immunology (IZI).

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.

# Citation

[![DOI](https://zenodo.org/badge/315346889.svg)](https://zenodo.org/badge/latestdoi/315346889)

BibTeX:
```
@software{dominik_j_otto_2020_4287137,
  author       = {Dominik J. Otto},
  title        = {fraunhofer-izi/tumor\_deconvolution: v1.0},
  month        = nov,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {v1.0},
  doi          = {10.5281/zenodo.4287137},
  url          = {https://doi.org/10.5281/zenodo.4287137} Titel anhand dieser DOI in Citavi-Projekt übernehmen
}
```
