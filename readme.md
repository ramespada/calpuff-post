# Post-processor for CALPUFF output files

> A python library for procesing CALPUFF (version > `7.2.1`) output files.


This library is capable of:
 - [ ] Reading CALPUFF output data files (`conc.dat`, `dryflx.dat`, `wetflx.dat`, `vis.dat`, etc.)
 - [ ] Get variable field for a given specie. 
 - [ ] Compute time-averaged magnitudes for a given specie and time-interval.
 - [ ] Get n-th higher values for a given specie and time-interval.

What this library does not do:
 - [ ] Make plots. (See `./scripts/plot_map.py`)
 - [ ] Make animations. (See `./scripts/make_anmation.py`)
 - [ ] Make tables. (See `./scripts/make_tables.py`)

## Dependencies:
This code uses the following python libraries:
 - `NumPy`
 - `Matplotlib` (optional, only for plots)
 - `Scipy` (optional, only for interpolate discrete values to regular grid)

## Usage:
```python

import calpost


puff = calpost.read_file('./data/conc_ducson.dat')

#Print file information:
puff.info()

# Get X,Y coordinates of receptors
X,Y = puff.get_coordinates()

# Get concentration of a specific pollutant
specie=puff.species[0]

C = puff.get_data(specie)

#Get time-averaged max value of a specific pollutant:

C_1hr = puff.get_data(specie, hours=1)

```

### Objects:

 - `COUT`: *Calpuff OUTput*
	+ `otype` (1: conc.dat, 2:dflx, 3:wflx, 4:dwflx, 5:vis)
 	+ `start_date`, `end_date`

### Functions:

 - `puff = calpost.read_file( <filepath>)`
 - `puff.info()`
 - `x, y = puff.get_coordinates(<specie>)`
 - `conc = puff.get_data(<specie>)`
 - `puff.get_time_avg_max(<specie>, <timestep_interval>)`


