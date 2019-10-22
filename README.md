# Python scripts to convert raw GPS data to estimated receiver height (snow depth) through GPS SNR

This is a work in progress.

# Required software:
  - https://www.unavco.org/software/data-processing/teqc/teqc.html  
  - https://github.com/kristinemlarson/gnssSNR

# Required data:
  - Raw GPS data (.lb2) - Your input
  - SP3 Orbital Data (.sp3) - ftp://cddis.nasa.gov/gnss/products/

Put the compiled or downloaded binaries in ./bin

# Steps
  - `process.py` for using tecq to extract RINEX .n and .o files from raw GPS data.
  - `snr.py` for using gnssSNR for extracting SNR from RINEX .o and SP3 file.
  - `analyze.py` for analysing the SNR data using elevation vs. SNR relationships through lomb scargle periodograms.
  - `plot.py` for plotting the analysed data.
