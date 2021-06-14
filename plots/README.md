 
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
[mass_limit_effect.py](mass_limit_effect.py): to see how the number of products blows up if the limit on the max MW allowed in the network is varied. The data used came from run logs, summarized below.

```bash
# This test was run on a machine with a Ryzen 5 4600H processor, 8 GB DDR4 RAM, GTX1650 GPU and 512GB SSD.
## MW limit =  200
Starting round 1
Took 0.06462359428405762 seconds to complete round 1
Products in this round: 16
Starting round 2
Took 0.7987122535705566 seconds to complete round 2
Products in this round: 110
Starting round 3
Took 26.78112268447876 seconds to complete round 3
Products in this round: 926
    ----
## MW Limit = 250
Starting round 1
Took 0.06683707237243652 seconds to complete round 1
Products in this round: 16
Starting round 2
Took 1.076136827468872 seconds to complete round 2
Products in this round: 148
Starting round 3
Took 117.04779601097107 seconds to complete round 3
Products in this round: 4537
    ----
## MW Limit = 300
Starting round 1
Took 0.06804275512695312 seconds to complete round 1
Original subset size: 16
Starting round 2
Took 1.0208632946014404 seconds to complete round 2
Original subset size: 184
Starting round 3
Took 571.1698200702667 seconds to complete round 3
Original subset size: 13832
```