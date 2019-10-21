# Python scripts for convert raw GPS data to RINEX (WIP)

This is a work in progress.

# Required libraries:
  - https://www.unavco.org/software/data-processing/teqc/teqc.html  
  - https://github.com/kristinemlarson/gnssSNR

Put the compiled or downloaded binaries in ./bin

# Steps
  - `process.py` for using tecq to extract RINEX .n and .o files from raw GPS data
  - `snr.py` for using gnssSNR for extracting SNR from RINEX .o and SP3 file
  - `analyze.py` for analysing the SNR data using elevation vs. SNR relationships through lomb scargle periodograms
