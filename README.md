# ReadASC

Read and plot ASC seismic ground motion records such as those obtained in the Turkey Earthquake.
Python code reads a set of 3 files one each for the NS, EW and Up or Z (vertical Component).
User interface uses tkinter.

Shows location of seismic instrument that originated the record on a map.
Plots acceleration, integrated velocity and displacement time history for each component.
Plot response spectra in SA vs Time Period or ADRS format - compute energy content of each component.
Plot 3D orbit plots.
Plot rotated resultant in the maximum acceleration, velocity or displacement directions.
Create RotD50, RotD00, RotD100 response spectra.  Compare to Geomean Spectra.
Plot resultant spectra in a Tripartite format.
Compare to ASCE 7-22 design spectra using any coordinates in the US - default is Downtown Los Angeles.
Save time vs. acceleration in a text format.

Example input files included. If desire is to read directly from the zip file containing all 3 components, use the branch ReadASCZip.

Needs many Python packages:
numpy
matplotlib
itertools
tkinter
scipy
