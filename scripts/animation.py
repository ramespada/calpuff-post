#!/usr/bin/python3

import math
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import animation
import calpost

def gen_levels(minval, maxval, nlevels):
    if minval <= 0 or maxval <= 0:
        raise ValueError("Both minval and maxval must be positive numbers")
    if minval >= maxval:
        raise ValueError("minval must be less than maxval")
    if nlevels < 1:
        raise ValueError("nlevels must be at least 1")

    # Calculate logarithmic range
    log_min = math.log10(minval)
    log_max = math.log10(maxval)

    # Generate initial candidates in log space
    log_candidates = []

    # Always include the endpoints
    log_candidates.append(log_min)
    log_candidates.append(log_max)

    # Generate key numbers (1 and 5 in each decade)
    current_decade = math.floor(log_min)
    while current_decade <= math.ceil(log_max):
        # Add 1 and 5 in this decade
        log1 = current_decade
        log5 = current_decade + math.log10(5)

        if log1 > log_min and log1 < log_max:
            log_candidates.append(log1)
        if log5 > log_min and log5 < log_max:
            log_candidates.append(log5)

        current_decade += 1

    # Add evenly spaced levels in log space (if needed to reach nlevels)
    if len(log_candidates) < nlevels + 1:
        additional_points = nlevels + 1 - len(log_candidates)
        step = (log_max - log_min) / (additional_points + 1)
        for i in range(1, additional_points + 1):
            log_point = log_min + i * step
            if log_point not in log_candidates:
                log_candidates.append(log_point)

    # Convert back to linear space and round to proper form
    candidates = [10**x for x in log_candidates]

    # Round to have only one non-zero digit
    def round_to_single_digit(x):
        if x == 0:
            return 0
        magnitude = 10 ** math.floor(math.log10(x))
        first_digit = round(x / magnitude)
        return first_digit * magnitude

    rounded = [round_to_single_digit(x) for x in candidates]

    # Remove duplicates and sort
    unique_levels = sorted(list(set(rounded)))

    # Ensure minval and maxval are exactly represented
    unique_levels[0] = minval
    unique_levels[-1] = maxval

    # If we still don't have enough levels, add more while maintaining constraints
    while len(unique_levels) < nlevels + 1:
        # Find the largest gap
        max_gap = 0
        gap_index = 0
        for i in range(len(unique_levels) - 1):
            gap = unique_levels[i+1] - unique_levels[i]
            if gap > max_gap:
                max_gap = gap
                gap_index = i

        # Add midpoint in log space
        log_mid = (math.log10(unique_levels[gap_index]) +
                   math.log10(unique_levels[gap_index+1])) / 2
        new_val = round_to_single_digit(10**log_mid)
        unique_levels.append(new_val)
        unique_levels = sorted(list(set(unique_levels)))

    return sorted(unique_levels)

# File path to calpuff output file
#file='./data/conc_mumbai.dat'
#file='./data/conc_losangeles.dat'
file='./data/conc_ducson.dat'
#file='./data/conc_waterloo.dat'

# Open and parse a CALPUFF output file
dat = calpost.read_file(file)

# Print access parsed information
dat.info()

# Get specie name
pollut=dat.species[0]

# Get discrete receptor concentrations.
X,Y   = dat.get_coordinates()

## Get gridded concentrations.
C     = dat.get_data(pollut)

#If not gridded receptors, then interpolate to regular grid.
if ( dat.ndrec > dat.ngrec ): # if discrete receptor is
    # Create a regular grid covering the domain of X and Y
    xi = np.arange(dat.x0, np.max(X), dat.dx)
    yi = np.arange(dat.y0, np.max(Y), dat.dy)
    X_grid, Y_grid = np.meshgrid(xi, yi)
    C_grid=np.zeros([dat.run_length, len(xi), len(yi)])
    points = np.column_stack((X, Y))

    # Perform bilinear (linear) interpolation
    for t in range(dat.run_length):
        C_grid[t,:,:] = griddata(points, C[t,:], (X_grid, Y_grid), method='linear')

    X=X_grid; Y=Y_grid
    C=C_grid*1e6

nlev=10
minval=max(1e-12, np.nanmin(C))
maxval=np.nanmax(C)

levels=gen_levels(minval,maxval,nlev)

fig, ax = plt.subplots()
contourf = [ax.contourf(X, Y, C[0,:,:], levels=levels, cmap='viridis',norm = LogNorm())]
cbar = plt.colorbar(contourf[0], ax=ax)

def update(t):
    for coll in contourf[0].collections:
        coll.remove()  # Remove old contours
    contourf[0] = ax.contourf(X, Y, C[t,:,:], levels=levels, cmap='viridis',norm = LogNorm())
    ax.set_title(f"Frame {t}")
    return contourf[0].collections

ani = animation.FuncAnimation(
    fig, update, frames=C.shape[0], blit=False
)

# Save as GIF using Pillow
ani.save("animation_"+pollut+".gif", writer='pillow', fps=5)

plt.close()
print("GIF saved as animation.gif")

