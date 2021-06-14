 
# Plot Scripts
This folder contains scripts that were used to create the plots used in our publications.

## FT-ICR-MS analysis
The code is contained in [mass_spectra_plot.ipynb](mass_spectra_plot.ipynb). A few notes:
- The raw data *.csv* files contained variable number of columns in some subsequent rows, which were causing trouble processing the data using pandas. As a quick resolution, we extracted columns of our interest and created new csv's using them. The original data can be found in the [data](../data/) directory.
- A sizable fraction of the data entries had not been assigned a molecular formula by the mass spectrometer (labelled "No Hit" in the "Molecular Formula" column). These peaks could be inorganic salt clusters. We cleaned the data by ignoring these peaks.
- The details of sample preparation, etc. have been covered in our manuscripts.

## Structural isomers counts
Count of structural isomers in the network as a function of generation, plotted against exact weight. Code is available in [struct_isomers.py](struct_isomers.py).

## Effect of upper limit on allowed molecular weight
[mass_limit_effect.py](mass_limit_effect.py): to see how the number of products blows up if the limit on the max MW allowed in the network is varied.
