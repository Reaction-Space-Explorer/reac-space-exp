 
# Plot Scripts
This folder contains scripts that were used to create the plots used in our publications.
## FT-ICR-MS analysis
The code is contained in [mass_spectra_plot.py](mass_spectra_plot.py). A few notes:
- The raw data *.csv* files contained variable number of columns in some subsequent rows, which were causing trouble processing the data using pandas.As a quick resolution, we extracted columns of our interest and created new csv's using them. The original data can be found in the [data](../data/) directory.