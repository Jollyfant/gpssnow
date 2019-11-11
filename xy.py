"""
Example script to plot satellite
"""

import os
import datetime
import matplotlib.pyplot as plt
import numpy as np

# Get collection of the SNR files
files = sorted(os.listdir("snr"))

from analyze import parseSNRFile, getPlaneCoordinates, extractMetadata

for file in files:

  filepath, date = extractMetadata(file)

  # Open a single data file for reading
  with open(filepath, "r") as infile:
    lines = infile.read().split("\n")[:-1]

  # Get all data points from the file
  data = list(map(lambda x: parseSNRFile(x, date), lines))

  for i in range(1, 33):

    if i != 20:
      continue

    # Filter currently active satellite
    SNRData = list(filter(lambda x: x["satellite"] == i, data))

    azi = np.array(list(map(lambda x: x["azimuth"], SNRData)))
    eli = np.array(list(map(lambda x: x["elevation"], SNRData)))

    # Get and show the coordinates
    x, y = getPlaneCoordinates(azi, eli)

    plt.scatter(x, y)
    plt.axhline(y=0)
    plt.axvline(x=0)
    plt.show()
