# HAMSTER-v2.0
Parallel version of the Heat And MoiSture Tracking framEwoRk

----

Updated version of the *Heat And MoiSture Tracking framEwoRk* (https://github.com/h-cel/hamster). It follows the WaterSip moisture tracking algorithm (https://gfi.wiki.uib.no/WaterSip) and adds bias correction with observed evaporation and precipitation, and tracking of sensible heat. The code has a new structure: two scripts (*preprocessing.py*, *tracking.py*), one with helper functions (*utilities.py*) and a control file (*namelist.input*). The tracking is parallelized with MPI.

## Repository layout

| Folder | Holds |
|---|---|
| `python/` | The reference code: `preprocessing.py`, `tracking.py`, `utilities.py`, `namelist.input`, `example/` (job scripts, sink mask) and `parallel-HAMSTER_forward/` (forward version). Copy the files to the run folder |
| `julia/` | `HAMSTER.jl`, a Julia port of preprocessing and tracking (backward and forward, MPI). Reads the raw data without overwriting it. See `julia/README.md` |

## Context

HAMSTER attributes the moisture sources of a precipitation event in a region (the sink region) by tracking the moisture gains and losses of the air parcels that contribute to it in the days before the event. Parcels follow the flow (Lagrangian view), so HAMSTER needs the output of a Lagrangian transport model. Only FLEXPART v11 output (https://www.flexpart.eu) is supported. Global FLEXPART-ERA5 trajectories for 1979-2024 with about 20 million parcels are available at https://zenodo.org/records/17952362. For another model, the input reading has to be adapted.

HAMSTER produces a two-dimensional field of moisture sources, like the one below:

<div align="center">
  <img src="https://github.com/user-attachments/assets/c59e8a29-7542-4607-a799-2d4f439220e6" width="500"/>
</div>

This field, E2P, gives for each grid cell the moisture evaporated there that ends up in the precipitation event (here the floods in Germany and Belgium in July 2021, black box). HAMSTER also produces the contribution of evaporation to integrated water vapor (E2Q) and of sensible heat to the temperature in the sink region (H2T).

## Components

| Script | Does |
|---|---|
| `preprocessing.py` | Reads the raw partoutputs and, for each pair of dates, computes the change in specific humidity (*Delta_q = qv[1,:] - qv[0,:]*) or potential temperature (*Delta_theta*) of every parcel. It corrects these changes with observed evaporation and precipitation: for example, a positive *Delta_q* in a cell with zero observed evaporation is set to 0, and a negative *Delta_q* without rain is set to 0. The correction acts en route, so it changes the relative contributions of the sources, not only the total (as in the previous version). The corrected values are written back to the partoutputs. Runs serially and only once |
| `tracking.py` | Reads the corrected partoutputs and the sink region, selects the parcels in the sink and applies WaterSip to get E2P, E2Q or H2T. It is run for each sink region. It tracks 30 days back: with `dt` = 3 h and 20 million parcels, *Delta_q* alone is 240 x 20,000,000, so it needs a lot of RAM and uses MPI |

Run `preprocessing.py` first, then `tracking.py`.

## Prerequisites

Python 3 with xarray, numba, netcdf4-python, mpi4py and h5py. On HPC-UGent:

```
module purge
module load xarray/2023.9.0-gfbf-2023a
module load numba/0.58.1-foss-2023a
module load netcdf4-python/1.6.4-foss-2023a
module load mpi4py/3.1.4-gompi-2023a
module load h5py/3.9.0-foss-2023a
```

## Running HAMSTER v2

Example submission scripts are in `python/example`. Preprocessing (serial):

```
#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l walltime=08:00:00
#PBS -l mem=6GB
#PBS -A 2022_204

python $run_path/preprocessing.py
```

Tracking (MPI):

```
#!/bin/bash

#PBS -l nodes=1:ppn=32
#PBS -l walltime=06:00:00
#PBS -l mem=120GB
#PBS -A 2022_204

mpirun -np $PBS_NP python $run_path/tracking.py
```

## Control file

```
[DATES]
start_date=set_dateini
end_date=set_datefin

[CONSTANTS]
num_days=30
dt=3
dx=0.5
dy=0.5
num_parc=19866694
atm_mass=set_mass

[PATHS]
path_data1=/XX/XX/set_path_to_partoutputs_folder
path_data2=/XX/XX/set_path_to_tracking_output
path_mask=/XX/XX/set_path_to_sink
path_obs=/XX/XX/set_path_to_observations

[NAMES]
obs_evap=evap_era5_3hourly
obs_heat=heat_era5_3hourly
obs_prec=prec_era5_3hourly

[FLAGS]
track_heat=false
obs_era5=true
save_e2q=true
save_diagnose=true
```

| Key | Meaning |
|---|---|
| `start_date`, `end_date` | Period of interest |
| `num_days` | Days to track the parcels |
| `dt`, `dx`, `dy` | Temporal (h) and spatial (degrees) resolution |
| `num_parc`, `atm_mass` | Number of parcels and total atmospheric mass in the simulation |
| `path_data1` | FLEXPART output (partoutputs) |
| `path_data2` | HAMSTER output |
| `path_mask` | Sink (or source) region file |
| `path_obs` | Observations for the bias correction |
| `obs_evap`, `obs_heat`, `obs_prec` | NetCDF file names of the observed evaporation, heat and precipitation |
| `track_heat` | Track heat |
| `obs_era5` | Observations are ERA5 |
| `save_e2q` | Save E2Q |
| `save_diagnose` | Save the E and P diagnostics |

## Backward and forward versions

The description above is the backward version: the sink region is given and the parcels that reach it during the rain event are tracked back. In the forward version the source region is given and the moisture that evaporates from it is tracked over the next 30 days. Its output fields are PFE (precipitation from evaporation), QFE (integrated water vapor from evaporation) and HFS (sensible heat from source). Example: PFE for North America, July 2021:

<div align="center">
  <img src="https://github.com/user-attachments/assets/020671cf-585c-4c89-b2fa-f64939225a23" width="500"/>
</div>

The forward code is in `python/parallel-HAMSTER_forward`.

## Caveats

- Only FLEXPART v11 output is read.
- Sink and source regions need the right coordinates and a variable named *mask*. Examples: `python/example` (sink) and `python/parallel-HAMSTER_forward/example` (source).
- Only ERA5 observations are supported. The observation time step must equal `dt` (3-hourly for `dt` = 3), and the dates must be centered (evaporation from 00 to 03 UTC is dated 01:30 UTC). Use for instance the CDO `settaxis` function, or change `preprocessing.py`.
- The tracking output file name has the date of the end of the period, but the time stamp inside is centered: the sources of rainfall between 00 and 03 UTC on 1 June 2015 are in `attribution_e2p_20150601030000.nc` with time 01:30 UTC.
- Know the code, not only how to run it, to avoid errors that go unnoticed.

## Referencing
If you use HAMSTER v2, please cite:

***Insua-Costa, D., Deman, V. M. H., and Miralles, D. G.*** Heat And MoiSture Tracking framEwoRk version 2.0 (HAMSTER v2) (2.0.0), Zenodo, https://doi.org/10.5281/zenodo.20181215, 2026.

***Deman, V. M. H., Insua-Costa, D., Keune, J., Koppa, A., and Miralles, D. G.*** A four-decade global Lagrangian air-parcel trajectory dataset for atmospheric moisture and heat analysis. Earth Syst. Sci. Data, 18, 4593–4615 (2026).

## References

***Sodemann, H., Schwierz, C., & Wernli, H.*** Interannual variability of Greenland winter precipitation sources: Lagrangian moisture diagnostic and North Atlantic Oscillation influence. J. Geophys. Res. D: Atmos., 113 (2008).

***Keune, J., Schumacher, D. L. & Miralles, D. G.*** A unified framework to estimate the origins of atmospheric moisture and heat using Lagrangian models. Geosci. Model Dev. 15, 1875–1898 (2022).

***Bakels, L., Tatsii, D., Tipka, A., Thompson, R., Dütsch, M., Blaschek, M., ... & Stohl, A.*** FLEXPART version 11: Improved accuracy, efficiency, and flexibility. Geosci. Model Dev., 17, 7595-7627 (2024).
