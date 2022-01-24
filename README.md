# MIMOSA
[![DOI](https://zenodo.org/badge/440905926.svg)](https://zenodo.org/badge/latestdoi/440905926)

MIMO Simulator with Amplifiers.
This repo simulates a Quadriga channel with GMP-based PAs and GMP-based DPD per antenna.

## Installation.
1. Download the repo.
2. Download Quadriga and extract the "quadriga_src" to the directory 'MIMOSA\src\modules\Channels'; .https://quadriga-channel-model.de/#Download

## Running.
The main script is in src/main.m.
The parameters of the test are in src/params.m. Use this file to adjust the PA model, etc.

## Citing.
If you use any part of this project, please cite it as:
```
@misc{TarverILADPD,
  author       = {Tarver, Chance},
  title        = {MIMOSA: MIMO Simulator with Amplifiers},
  month        = jan,
  year         = 2022,
  doi          = {put apprpriate doi here from current doi above},
}
```
