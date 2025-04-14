# Post-processor for CALPUFF output files

> A python library for procesing CALPUFF (version > `7.2.1`) output files.

This library is capable of:
 - [x] Reading CALPUFF output files (`conc.dat`, `dryflx.dat`, `wetflx.dat`, `vis.dat`, etc.)
 - [x] Get output variable values for a given specie. 
 - [X] Get n-th higher values for a given specie and time-interval.

What this library does NOT do:
 - Make plots. (See `./tests/plot_avg.py`)
 - Make animations. (See `./tests/animation.py`)
 - Make tables. (See `./tests/tables.py`)

## Dependencies:
This code uses the following python libraries:
 - `NumPy`

## Usage:
```python

import calpost

puff = calpost.read_file('./data/conc.dat')

#Print file information:
puff.info()

# Get X,Y coordinates of receptors
X,Y = puff.get_coordinates()

# Get concentration of a specific pollutant

specie = puff.species[0]
C = puff.get_data(specie)

# Get time-averaged max value of a specific pollutant:

C_1hr = puff.get_time_avg_max(specie, interval=1, rank=1)
```

### Objects:
 - `CalpuffOutput`

### Functions:

 - `puff = calpost.read_file( <filepath>)`: read file and stores in `puff` variables that describes the dataset.
 - `puff.info()`: prints header file information.
 - `x, y = puff.get_coordinates(<specie>)`: return 2 arrays with the x,y coordinates velues for each receptor.
 - `conc = puff.get_data(<specie>)`: returns an array with the main magnitude store in the file (concentration, fluxes, etc.)
 - `puff.get_time_avg_max(<specie>, <timestep_interval>, <rank>)`: returns the averaged nth-max values of the main magnitude for a given specie, interval and rank.

### To-do list:
- [ ] Implement the rank tables function that returns the nth maximum values found.
- [ ] add an optional argument "group" for `get_data` and `get_coordinates`.
- [ ] test it with other output files types (dflx.dat, wflx.dat, vis.dat, etc.)
- [ ] test it with complex terrain receptors.
- [ ] purpose new file structure for newer calpuff versions to come.
- [ ] document, and add descriptions to all functions with the expected arguments and output types.
