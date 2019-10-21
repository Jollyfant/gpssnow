"""

Python script to infer snow height from GPS reflective interference
Modified after MatLab script "analyse.m" provided by Flavio Canavo (INGV)

Author: Mathijs Koymans, 2019, KNMI

"""

# Load the required libraries
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

from dateutil.parser import parse
from scipy.signal import lombscargle

def parseSNRFile(line):

  """
  Def parseSNRFile
  Reads out an SNR file created by program: gnssSNR
  """

  # Columns are defined in: https://github.com/kristinemlarson/gnssSNR
  parameters = line.split()

  return {
    "satellite": int(parameters[0]),
    "elevation": float(parameters[1]),
    "azimuth": float(parameters[2]),
    "seconds": float(parameters[3]),
    "reflector": float(parameters[4]),
    "S6": float(parameters[5]),
    "S1": float(parameters[6]),
    "S2": float(parameters[7]),
    "S5": float(parameters[8]),
    "S7": float(parameters[9]),
    "S8": float(parameters[10])
  }


def createPlot(x, y, s, z, f, p):

  """
  Def createPlot
  Creates a plot with three subplots
   [1] Original SNR data with 2nd order polynomial fitted (orange)
   [2] Fitted polynomial subtracted from SNR data
   [3] Lomb-Scargle periodogram of the remaining data
  """

  # One
  plt.subplot(3, 1, 1)
  plt.scatter(x, s, s=1)
  plt.plot(x, z, c="orange")

  # Two
  plt.subplot(3, 1, 2)
  plt.scatter(x, y, s=1)

  # Three
  plt.subplot(3, 1, 3)
  plt.plot(f, p)

  # Plot the maximum
  ind = np.argmax(p)

  maxP = p[ind]
  maxF = f[ind]

  if np.mean(p) + 3 * np.std(p) < maxP:
    color = "green"
  else:
    color = "red"

  plt.scatter(maxF, maxP, c=color)

  plt.show()


def getPlaneCoordinates(azimuth, elevation):

  """
  def getPlaneCoordinates
  Returns x, y plane coordinates based on azimuth and elevation
  angles with respect to the GPS receiver
  """

  # In meters
  RECEIVER_HEIGHT = 2

  # Height of the receiver divided by tangent of angle equals
  # the distance of the receiver
  r = RECEIVER_HEIGHT * np.tan(elevation)

  # Combine with the azimuth to find the x, y coordinate of the reflection
  return r * np.sin(azimuth), r * np.cos(azimuth) 


def extractMetadata(file):

  """
  def extractMetadata
  Extracts metadata from a SNR file
  """

  # Parse the filename to get info
  return os.path.join("snr", file), parse(file[4:-10])


def spatialFilter(obj):

  """
  def spatialFilter
  May filter reflections based on their spatial coordinates
  """

  return True
  # Calculate the relative x, y coordinates
  x, y = getPlaneCoordinates(obj["azimuth"], obj["elevation"])

  return x >= 0 and y >= 0

def processSNRFile(data, date):

  """
  def processSNRFile
  Processes a single SNR file to extract snow height information
  """

  NUMBER_OF_SATELLITES = 32
  NUMBER_OF_SAMPLES_MIN = 2000

  # Go over all possible satellites in the file (32)
  for i in range(NUMBER_OF_SATELLITES):
  
    print(i)
    # Filter currently active satellite
    datas = filter(lambda x: x["satellite"] == i, data)
  
    # Implement spatial filter based on receiver height
    datas = list(filter(spatialFilter, datas))
  
    # No data after filters
    if len(datas) == 0:
      continue
  
    # Get the S1 and Elevation parameters
    E = np.array(list(map(lambda x: x["elevation"], datas)))
    S1 = np.array(list(map(lambda x: x["S1"], datas)))
  
    # Stop when elevation begins to decrease and satellite has passed overhead
    # Can we used the removed information?
    ind = np.argmax(np.diff(E) <= 0)
    E = E[:ind]
    S1 = S1[:ind]
  
    # Sort by elevation: extract indices
    idx = E.argsort()
    E = E[idx]
    S1 = S1[idx]
  
    # Skip anything with not enough samples
    if len(E) < NUMBER_OF_SAMPLES_MIN:
      continue
  
    # Fit a 2nd order polynomial to the data
    polynomial = np.poly1d(np.polyfit(E, S1, 2))
  
    # Apply and subtract the fitted polynomial
    dS1 = S1 - polynomial(E)
  
    # Convert decibels to volts?
    vdSNR = (10.0 ** (dS1 / 10.0))
  
    # Lomb scarlge algorithm to get periods
    dE = np.sin(np.radians(E))
  
    # List of frequencies to get the power at
    freqs = np.linspace(0.001, 1000, 1E3)
  
    # Get the power at different periods
    power = lombscargle(dE, dS1, freqs)
  
    # Create the plot
    createPlot(E, dS1, S1, polynomial(E), freqs, power)


if __name__ == "__main__":

  # Go over the directory with gnssSNR files
  for file in sorted(os.listdir("snr")):

    # Create the filepath and extract the file date from its filename
    filepath, date = extractMetadata(file)

    # Open a single data file for reading
    with open(filepath, "r") as infile:
      lines = infile.read().split("\n")[:-1]

    data = list(map(parseSNRFile, lines))

    processSNRFile(data, date)
