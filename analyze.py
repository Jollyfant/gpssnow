"""

Python script to infer snow height from GPS reflective interference
Modified after MatLab script "analyse.m" provided by Flavio Canavo (INGV)

See Larson & Nievinski, 2012

Author: Mathijs Koymans, Oct 2019, KNMI

"""

# Load the required libraries
import random
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

from dateutil.parser import parse
from scipy.signal import lombscargle

SHOW_PLOT = True
SPATIAL_FILTER = False

def parseSNRFile(line):

  """
  Def parseSNRFile
  Reads out an SNR file created by program: gnssSNR
  """

  # Columns are defined in: https://github.com/kristinemlarson/gnssSNR
  parameters = line.split()

  # Important parameters are elevation, azimuth, S1 (L1)
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


def freq2height(freq):

  """
  def freq2height
  Converts peak SNR frequency to receiver height
  """

  # GPS frequencies in Hz
  L1_GPS_FREQUENCY = 1575.42E6
  L2_GPS_FREQUENCY = 1227.60E6

  # Wavelength = c / f
  LAMBDA = 299792458. / L1_GPS_FREQUENCY

  return (LAMBDA * freq) / (4. * np.pi)


def getHeight(freqs, power):

  """
  def getHeight
  Returns the frequency of the peak of the lomb scargle periodogram
  """

  # Get the index of maximum power
  ind = np.argmax(power)

  # Extract maximum power and frequency
  pMax = power[ind]
  fMax = freqs[ind]

  # Extract h from equation (2)
  return freq2height(fMax)


def createPlot(x, y, s, z, f, p, date, i):

  """
  Def createPlot
  Creates a plot with three subplots
   [1] Original SNR data with 2nd order polynomial fitted (orange)
   [2] Fitted polynomial subtracted from SNR data
   [3] Lomb-Scargle periodogram of the remaining data
  """

  # Add date to the plot
  plt.suptitle(date.isoformat() + " - Satellite: " + str(i))

  # One
  ax = plt.subplot(3, 1, 1)
  ax.set_ylabel("Volts")
  plt.scatter(x, s, s=1)
  plt.plot(x, z, c="orange")

  # Two
  ax = plt.subplot(3, 1, 2)
  ax.set_ylabel("Volts")
  ax.set_xlabel("Elevation Angle")
  plt.scatter(x, y, s=1)

  # Three
  ax = plt.subplot(3, 1, 3)
  ax.set_xlabel("Reflector Height (m)")
  ax.set_ylabel("Relative Power")

  # Get the maximum power pxx
  ind = np.argmax(p)

  pMax = p[ind]
  fMax = f[ind]

  # Plot reflector height (m) instead of angular frequency (Hz)
  f = freq2height(f)
  fMax = freq2height(fMax)

  # Plot the lomb scargle spectrum
  plt.plot(f, p)

  # Plot the peak value (green if significant)
  if np.mean(p) + 2 * np.std(p) < pMax:
    color = "green"
  else:
    color = "red"

  # Show the peak
  plt.scatter(fMax, pMax, c=color)

  plt.show()


def getPlaneCoordinates(azimuth, elevationAngle):

  """
  def getPlaneCoordinates
  Returns x, y plane coordinates based on azimuth and elevation
  angles with respect to the height of the GPS receiver
  """

  # In meters
  RECEIVER_HEIGHT = 2

  # Height of the receiver divided by tangent of angle equals
  # the distance of the receiver
  r = RECEIVER_HEIGHT * np.tan(elevationAngle)

  # Combine with the azimuth to find the x, y coordinate of the reflection
  return r * np.sin(azimuth), r * np.cos(azimuth) 


def extractMetadata(file):

  """
  def extractMetadata
  Extracts metadata from a SNR file
  """

  # Parse the filename to get info
  return (
    os.path.join("snr", file),
    parse(file[4:-10])
  )


def spatialFilter(obj):

  """
  def spatialFilter
  May filter reflections based on their spatial coordinates
  """

  # Calculate the relative x, y coordinates
  x, y = getPlaneCoordinates(obj["azimuth"], obj["elevation"])

  return x >= 0 and y >= 0

def sanitize(E, S1):

  """
  def sanitize
  Sanitizes the input data and prepares it for processing
  """

  # Stop when elevation begins to decrease and satellite has passed overhead
  # Can we used the removed information?
  ind = np.argmax(np.diff(E) <= 0)
  E = E[:ind]
  S1 = S1[:ind]
 
  # Sort by elevation: extract indices
  idx = E.argsort()
  E = E[idx]
  S1 = S1[idx]

  return E, S1

def extract(SNRData):

  """
  def extract
  Extracts the elevation, S1 data from a parsed SNR file
  """
  # Get the S1 and Elevation parameters

  E = np.array(list(map(lambda x: x["elevation"], SNRData)))
  S1 = np.array(list(map(lambda x: x["S1"], SNRData)))

  return E, S1

def processSNRFile(data, date):

  """
  def processSNRFile
  Processes a single SNR file to extract snow height information
  """

  NUMBER_OF_SATELLITES = 32
  NUMBER_OF_SAMPLES_MIN = 2000

  # Go over all possible satellites in the file (32)
  for i in range(NUMBER_OF_SATELLITES):
  
    # Filter currently active satellite
    SNRData = list(filter(lambda x: x["satellite"] == i, data))
  
    # Implement spatial filter based on receiver height
    if SPATIAL_FILTER:
      SNRData = list(filter(spatialFilter, SNRData))
  
    # No data after filters
    if len(SNRData) == 0:
      continue
  
    # Extract the elevation and S1 (L1) data
    E, S1 = extract(SNRData)
  
    # Clean the data
    E, S1 = sanitize(E, S1)
  
    # Skip anything with not enough samples
    if len(E) < NUMBER_OF_SAMPLES_MIN:
      continue
  
    # Convert dBV to Volts: unclear to what reference the dB is expressed
    # According to Shuanggen et al., 2016 this is correct:
    vS1 = 10.0 ** (S1 / 20.0)

    # Fit a 2nd order polynomial to the voltage data
    # Remove Pd & Pr (direct signal)
    polynomial = np.poly1d(np.polyfit(E, vS1, 2))
  
    # Apply and subtract the fitted polynomial
    dS1 = vS1 - polynomial(E)
  
    # From Larson & Nievinski, 2012 - equation (2)
    # 
    # SNR ~ A * cos(4 * np.pi * h * (lambda ** -1) * sin(e) + phi) (eq. 2)
    # 
    # To find h, we have to identify the frequency of a function SNR = a * cos(bx + c)
    #
    # Where:
    #    b = 4 * np.pi * h * (lambda ** -1)
    #    x = sin(e) (elevation)
    # 
    # The amplitude (A, a), phase shift (c, phi) can be ignored

    # Create sin(e) by taking the sine of the elevation
    sE = np.sin(np.radians(E))

    # List of angular frequencies to get the power at: should capture the peak
    # 
    freqs = np.linspace(0.01, 1000, 1E4)

    # Get the power at different periods using the Lomb Scargle algorithm for unregularly sampled data
    # Look at frequency content of SNR (dS1) as a function of sin(elevation)
    power = lombscargle(
      sE,
      dS1,
      freqs,
      normalize=True
    )
   
    # Create the plot
    if SHOW_PLOT:
      createPlot(E, dS1, vS1, polynomial(E), freqs, power, date, i)

    # Get the reflection height
    h = getHeight(freqs, power)

    # Write the date & expected receiver height to an outfile
    with open("outfile.csv", "a") as outfile:
      outfile.write(date.isoformat() + " " + str(h) + "\n")

if __name__ == "__main__":

  # Go over the directory with gnssSNR files
  for file in sorted(os.listdir("snr")): 

    # Create the filepath and extract the file date from its filename
    filepath, date = extractMetadata(file)

    # Open a single data file for reading
    with open(filepath, "r") as infile:
      lines = infile.read().split("\n")[:-1]

    # Get all data points from the file
    data = list(map(parseSNRFile, lines))

    processSNRFile(data, date)
