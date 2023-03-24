# MIMOSA
[![DOI](https://zenodo.org/badge/440905926.svg)](https://zenodo.org/badge/latestdoi/440905926)

MIMO Simulator with Amplifiers (MIMOSA).
This repo simulates a Quadriga channel with GMP-based PAs and GMP-based DPD per antenna.

## Installation
1. Download the repo.
2. Download Quadriga and extract the "quadriga_src" to the directory 'MIMOSA\src\modules\Channels'; .https://quadriga-channel-model.de/#Download

## Quick Start
Firstly add the "modules" directory to your matlab path. 
The main script is in applications/GMP_MIMO/main.m. The parameters of the test are in params.m. Use this file to adjust the PA model, etc.

## Publications

A version of this repo was used to create: 
C. Tarver, A. Balatsoukas-Stimming, C. Studer and J. R. Cavallaro, "[OFDM-Based Beam-Oriented Digital Predistortion for Massive MIMO](https://ieeexplore.ieee.org/document/9401479)," 2021 IEEE International Symposium on Circuits and Systems (ISCAS), Daegu, Korea, 2021, pp. 1-5, doi: 10.1109/ISCAS51556.2021.9401479.


## Using
The main idea behind this repo is to create an object-oriented library of DSP modules that are independednt of any specific project or application. Using these modules, we can create various applications to quickly prototype new ideas out in simulation. 

## Citing.
If you use any part of this project, please cite it as:
```
@misc{tarver_mimosa,
  author       = {Tarver, Chance},
  title        = {MIMOSA: MIMO Simulator with Amplifiers},
  month        = jan,
  year         = 2022,
  doi          = {put apprpriate doi here from current doi above},
}
```
