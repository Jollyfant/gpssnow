"""

Python script to infer snow height from GPS reflective interference
Modified after MatLab script "analyse.m" provided by Flavio Canavo (INGV)

See Larson & Nievinski, 2012

Author: Mathijs Koymans, Oct 2019, KNMI

"""

# Load the required libraries
import numpy as np
import os
import signal
import matplotlib.pyplot as plt
import datetime
import multiprocessing 

from dateutil.parser import parse
from scipy.signal import lombscargle

# Some options
SHOW_PLOT = False
SPATIAL_FILTER = True
APPLY_CORR = True
POLY_ORDER = 2

def parseSNRFile(line, date):

  """
  Def parseSNRFile
  Reads out an SNR file created by program: gnssSNR
  """

  # Columns are defined in: https://github.com/kristinemlarson/gnssSNR
  parameters = line.split()

  elevation = float(parameters[1])

  # Apply P/T correction (Peng et al., 2019)
  # Application of GNSS interferometric reflectometry for detecting storm surge
  if APPLY_CORR:
    elevation += getCorrection(elevation, date)

  # Important parameters are elevation, azimuth, S1 (L1)
  return {
    "satellite": int(parameters[0]),
    "elevation": elevation,
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

def freq2height(freq, type):

  """
  def freq2height
  Converts peak SNR frequency to receiver height
  """

  # GPS frequencies in Hz
  L1_GPS_FREQUENCY = 1575.42E6
  L2_GPS_FREQUENCY = 1227.60E6

  # Wavelength = c / f
  LAMBDA1 = 299792458. / L1_GPS_FREQUENCY
  LAMBDA2 = 299792458. / L2_GPS_FREQUENCY

  if type == "S1":
    return (LAMBDA1 * freq) / (4. * np.pi)
  elif type == "S2":
    return (LAMBDA2 * freq) / (4. * np.pi)

def getTemperature(date):

  """
  def getTemperature
  Returns the temperature that belongs to a given date (month)
  """

  # Average monthly temperatures
  MONTHLY_TEMPERATURES = [
    1, -1, 5, 11, 18, 22, 23, 21, 17, 12, 10, 8
  ]

  # Month index starts at 1
  return MONTHLY_TEMPERATURES[date.month - 1]

def getCorrection(elevations, date):

  """
  def getCorrection
  Elevation height (degrees) - The Calculation of Astronomical Refraction in Marine Navigation
  DOI: https://doi.org/10.1017/S0373463300022037
  TODO:: implement H
  """

  # hPa at 1800 meters altitude
  PRESSURE = 800
  temperature = getTemperature(date) 

  phi = np.radians(elevations + (7.31 / (elevations + 4.4)))

  return (510. / ((9 * temperature / 5) + 492)) * (PRESSURE / 1010.16) * (1. / np.tan(phi))

def getHeight(freqs, power, type):

  """
  def getHeight
  Returns the frequency of the peak of the lomb scargle periodogram
  """

  # Get the index of maximum power
  ind = np.argmax(power)

  # Extract maximum power and frequency
  fMax = freqs[ind]

  # Extract h from equation (2)
  return freq2height(fMax, type)


def createPlot(x, y, s, z, f, p, date, i, type):

  """
  Def createPlot
  Creates a plot with three subplots
   [1] Original SNR data with 2nd order polynomial fitted (orange)
   [2] Fitted polynomial subtracted from SNR data
   [3] Lomb-Scargle periodogram of the remaining data
  """

  # Add date to the plot
  plt.suptitle(date.isoformat() + " - Satellite: " + str(i) + " " + type)

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
  f = freq2height(f, type)
  fMax = freq2height(fMax, type)

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
  plt.clf()


def getPlaneCoordinates(azimuth, elevationAngle):

  """
  def getPlaneCoordinates
  Returns x, y plane coordinates based on azimuth and elevation
  angles with respect to the height of the GPS receiver
  """

  azimuth = np.radians(azimuth)
  elevationAngle = np.radians(elevationAngle)

  # GPS receiver height in meters
  RECEIVER_HEIGHT = 1.5

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

  # Sort by elevation: extract indices
  idx = E.argsort()
  E = E[idx]
  S1 = S1[idx]

  return E, S1

def extract(SNRData, key, date):

  """
  def extract
  Extracts the elevation, S1 data from a parsed SNR file
  """

  # Get the S1 and Elevation parameters
  E = np.array(list(map(lambda x: x["elevation"], SNRData)))
  S1 = np.array(list(map(lambda x: x[key], SNRData)))
  s = np.array(list(map(lambda x: x["seconds"], SNRData)))

  return E, S1, s

def split(E, S1, s):

  traces = list()

  while True:

    ind = np.argmax(np.diff(s) > 3600)

    if ind == 0:
      break

    # Save one trace
    E1 = E[:ind]
    S11 = S1[:ind]
    s1 = s[:ind]

    traces.append((E1, S11, s1))

    # Cut the next
    E = E[ind + 1:]
    S1 = S1[ind + 1:]
    s = s[ind + 1:]

  # Remember the remaining trace!
  traces.append((E, S1, s))

  return traces
  
def validate(E, S1, s):

  """
  def validate
  Validates E, S1 by some simple quality metrics
  """

  # Minimum number of samples
  NUMBER_OF_SAMPLES_MIN = 2000

  # Get the length of the trace
  # Skip anything with not enough samples
  if len(E) < NUMBER_OF_SAMPLES_MIN:
    return False

  # Remove traces with high scatter
  if np.std(np.diff(S1)) > 1:
    return False

  # Some weird hyperbole?!
  _, r, _, _, _ = np.polyfit(s, E, 1, full=True)
  if r[0] / len(s) > 0.1:
    return False

  return True

def processSNRFile(data, date):

  """
  def processSNRFile
  Processes a single SNR file to extract snow height information
  """

  results = {
    "date": date,
    "values": list()
  }

  # Go over all possible GPS satellites in the snr file (1 - 32)
  for i in range(1, 33):
  
    # Filter currently active satellite
    SNRData = list(filter(lambda x: x["satellite"] == i, data))
  
    # Implement spatial filter based on receiver height
    if SPATIAL_FILTER:
      SNRData = list(filter(spatialFilter, SNRData))
  
    # No data after filters
    if len(SNRData) <= 1:
      continue

    # Extract the elevation and S1 (L1) data
    for (E, S1, s) in split(*extract(SNRData, "S1", date)):
      height = processSignal(i, E, S1, s, "S1", date)
      if height is not None:
        results["values"].append((height, i, "S1"))

    # Extract the elevation and S1 (L1) data
    for (E, S2, s) in split(*extract(SNRData, "S2", date)):
      height = processSignal(i, E, S2, s, "S2", date)
      if height is not None:
        results["values"].append((height, i, "S2"))

  return results

def processSignal(i, E, S1, s, type, date):

  # Clean the data
  #E, S1 = sanitize(E, S1)
 
  # Validate the trace for its quality
  if not validate(E, S1, s):
    return None

  # Convert dBV to Volts: unclear to what reference the dB is expressed
  # According to Shuanggen et al., 2016 this is correct:
  vS1 = 10.0 ** (S1 / 20.0)

  # Fit a 2nd order polynomial to the voltage data
  # Remove Pd & Pr (direct signal)
  polynomial = np.poly1d(np.polyfit(E, vS1, POLY_ORDER))

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

  # Linear space of angular frequencies to get the power at: should capture the peak
  # Between 
  freqs = np.linspace(1, 200, 1E4)

  # Get the power at different periods using the Lomb Scargle algorithm for unregularly sampled data
  # Look at frequency content of SNR (dS1) as a function of sin(elevation)
  power = lombscargle(
    sE,
    dS1,
    freqs
  )
  
  # Create the plot
  if SHOW_PLOT:
    createPlot(E, dS1, vS1, polynomial(E), freqs, power, date, i, type)

  # Get the reflection height
  return getHeight(freqs, power, type)

def worker(file):

  """
  def worker
  Worker function passed to multiprocessing unit
  """

  # Create the filepath and extract the file date from its filename
  filepath, date = extractMetadata(file)

  if date < datetime.datetime(2017, 3, 25):
    return

  # Open a single data file for reading
  with open(filepath, "r") as infile:
    lines = infile.read().split("\n")[:-1]

  # Get all data points from the file
  data = list(map(lambda x: parseSNRFile(x, date), lines))

  return processSNRFile(data, date)

def initWorker():

  """
  def initWorker
  Initialzes a worker
  """

  print("Initializing new worker process: %s." % multiprocessing.current_process())

  # Ignore signal interrupts in the worker process
  signal.signal(signal.SIGINT, signal.SIG_IGN)

if __name__ == "__main__":

  """
  def __main__
  Extracts reflector height using GNSS reflective interferometry  
  """

  print("Initializing pool of %i workers for processing." % multiprocessing.cpu_count())

  # Get collection of the SNR files
  files = sorted(os.listdir("snr"))

  if SHOW_PLOT:
    for file in files:
      worker(file)

  # Create a pool (one process for each core)
  pool = multiprocessing.Pool(multiprocessing.cpu_count(), initWorker)

  try:
    # Use a pool of four workers
    for result in pool.imap(worker, files):

      if result is None:
        continue

      print("Completed SNR process for %s." % result["date"].isoformat())

      # Write each result to a line
      for (height, i, type) in result["values"]:
        with open("outfile.dat", "a") as outfile:
          outfile.write("%s %.3f %i %s\n" % (result["date"].isoformat(), height, i, type))

  except KeyboardInterrupt:
    pool.terminate()

  finally:
    pool.close()
    pool.join()
