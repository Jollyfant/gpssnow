#!/usr/local/bin/python3

"""

Python script to infer snow height from GPS reflective interference
Modified after MatLab script "analyse.m" provided by Flavio Canavo (INGV)

Implements the approach of Larson & Nievinski, 2012

Author: Mathijs Koymans, Nov 2019, KNMI

"""

# Load the required libraries
import numpy as np
import os
import datetime
import multiprocessing
import signal

from dateutil.parser import parse

# Hoist configuration to global scope
configuration = None

def printVerbose(message):

  """
  def printVerbose
  Prints verbose message
  """
  if configuration.v:
    print(datetime.datetime.now().isoformat(), message)

def parseSNRFile(filepath, date):

  """
  def parseSNRFile
  Parses SNR File
  """

  # Read the file
  with open(filepath, "r") as infile:
    lines = infile.read().split("\n")[:-1]

  satellites = list()
  elevations = list()
  azimuths = list()
  seconds = list()
  S1s = list()
  S2s = list()

  # Parse the file contents
  for line in lines:

    # Each line
    (satellite, elevation, azimuth, second, _, _, S1, S2, *void) = line.split()

    satellites.append(int(satellite))
    elevations.append(float(elevation))
    azimuths.append(float(azimuth))
    seconds.append(float(second))
    S1s.append(float(S1))
    S2s.append(float(S2))

  # Apply the correction for atmosphere if requested
  if configuration.correction:
    elevations += getCorrection(np.array(elevations), date)

  # Create numpy arrays
  data = dict({
    "satellites": np.array(satellites),
    "elevations": np.array(elevations),
    "azimuths": np.array(azimuths),
    "seconds": np.array(seconds),
    "S1": np.array(S1s),
    "S2": np.array(S2s)
  })

  return data

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
  """

  # hPa at 1800 meters altitude
  #PRESSURE = 800
  # hPa at 2500 meters altitude
  PRESSURE = 709
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

  return freq2height(fMax, type)

def createPlot(x, y, s, z, f, p, type):

  """
  Def createPlot
  Creates a plot with three subplots
   [1] Original SNR data with 2nd order polynomial fitted (orange)
   [2] Fitted polynomial subtracted from SNR data
   [3] Lomb-Scargle periodogram of the remaining data

   if date != datetime.datetime(2015, 4, 3):
     return None
  """

  import matplotlib
  import matplotlib.ticker as plticker
  import matplotlib.pyplot as plt

  matplotlib.rcParams["mathtext.fontset"] = "stix"
  matplotlib.rcParams["font.family"] = "STIXGeneral"
  matplotlib.rcParams.update({"font.size": 12})

  x = np.sin(np.radians(x))

  # Get the maximum power pxx
  ind = np.argmax(p)

  pMax = p[ind]
  fMax = f[ind]

  # Plot reflector height (m) instead of angular frequency (Hz)
  h = freq2height(f, type)
  hMax = freq2height(fMax, type)

  # One
  ax = plt.subplot(2, 1, 1)
  ax.set_ylabel("SNR (Volts)")
  ax.set_ylim(50, 250)
  plt.plot(x, s, label="SNR data")
  plt.plot(x, z, "--", c="orange", label="Polynomial fit")
  plt.legend()

  ax.yaxis.set_minor_locator(plticker.MultipleLocator(base=10))
  ax.xaxis.set_major_locator(plticker.MultipleLocator(base=5))
  ax.xaxis.set_minor_locator(plticker.MultipleLocator(base=0.5))
  ax.xaxis.set_major_formatter(plt.NullFormatter())

  # Two
  ax = plt.subplot(2, 1, 2)
  ax.set_ylabel("SNR (Volts)")
  ax.set_xlabel("Sine Elevation Angle")
  plt.plot(x, y)
  plt.plot(x, 50 * np.cos(fMax * x + 1.2 * np.pi), "--", c="red", label="Dominant Frequency")
  ax.set_ylim(-100, 100)
  plt.legend()

  #ax.yaxis.set_minor_locator(plticker.MultipleLocator(base=10))
  #ax.xaxis.set_major_locator(plticker.MultipleLocator(base=5))
  #ax.xaxis.set_minor_locator(plticker.MultipleLocator(base=0.5))

  fig = plt.gcf()
  fig.set_size_inches(10.5, 0.5 * 10.5)
  fig.savefig("volts.pdf", dpi=100, bbox_inches="tight")

  plt.show()
  plt.clf()

  # Lomb Scargle periodogram
  ax = plt.subplot(1, 1, 1)
  ax.set_xlabel("Reflector height (m)")
  ax.set_ylabel("Normalized power")

  ax.set_ylim(0, 0.7)
  ax.set_xlim(0, 3)
  ax.yaxis.set_minor_locator(plticker.MultipleLocator(base=0.05))
  ax.xaxis.set_minor_locator(plticker.MultipleLocator(base=0.1))

  # Plot the lomb scargle spectrum
  plt.plot(h, p)

  # Plot the peak value (green if significant)
  if np.mean(p) + 2 * np.std(p) < pMax:
    color = "green"
  else:
    color = "red"

  # Show the peak
  plt.scatter(hMax, pMax, c=color, label="Dominant frequency $\omega_{\mathrm{P_{max}}}$")
  plt.legend()
  plt.show()

  fig = plt.gcf()
  fig.set_size_inches(10.5, 0.25 * 10.5)
  fig.savefig("periodogram.pdf", dpi=100, bbox_inches="tight")


def getPlaneCoordinates(elevationAngle, azimuth):

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
  return r * np.cos(azimuth), r * np.sin(azimuth) 

def spatialFilter(azimuths):

  """
  def spatialFilter
  May filter reflections based on their spatial coordinates
  """

  # Create a bitwise mask
  mask = np.ones(len(azimuths), dtype=bool)

  # Bitwise AND must be in between for all filters
  for values in configuration.azimuths:
    mask &= azimuths > values[0]
    mask &= azimuths < values[1]

  # Check indices that pass
  return np.where(mask)


def splitTraces(seconds):

  """
  def splitTraces
  Splits a satellite trace in multiple overhead passes
  Tried splitting by ascending/descending traces: end up with problems
  for satellites that have reflections curved in the xy-plane
  """

  traces = list()
  prev = 0

  while True:

    # A gap of over 1000s should be a new pass
    # Index of first gap exceeding this value
    ind = np.argmax(np.diff(seconds) > 3600)

    # Nothing found: stop!
    if ind == 0:
      break

    # Save this particular trace
    traces.append((prev, ind))
    prev = ind + 1

    # Cut the seconds
    seconds = seconds[prev:]

  # The remaining trace of course (with no end)
  traces.append((prev, None))

  return traces
  

def validateTrace(E):

  """
  def validateTrace
  Validates E, S1 by some simple quality metrics
  """

  # Get the length of the trace
  # Skip anything with not enough samples
  if len(E) < configuration.samples:
    return False

  # The satellite elevation angle must monotonically increase or decrease
  if not (np.all(np.diff(E) < 0) or np.all(np.diff(E) > 0)):
    return False

  return True


def processSNRFile(data):

  """
  def processSNRFile
  Processes a single SNR file to extract snow height information
  """

  # Create a container for the results
  results = list()

  # Go over all available GPS satellites
  for i in configuration.satellites:
  
    # Extract data for the correct satellite
    idx = data["satellites"] == i

    # Implement spatial filter based on receiver height
    elevations = data["elevations"][idx]
    azimuths = data["azimuths"][idx]
    seconds = data["seconds"][idx]
    S1 = data["S1"][idx]
    S2 = data["S2"][idx]
  
    # Implement spatial filter based on receiver height
    if configuration.azimuths is not None:

      # Make sure to implement this for your azimuths
      idx = spatialFilter(azimuths)

      elevations = elevations[idx]
      seconds = seconds[idx]
      azimuths = azimuths[idx]
      S1 = S1[idx]
      S2 = S2[idx]

    # No data after spatial filter: must have more than one sample for coming np.diff
    if len(seconds) <= 1:
      continue

    # Process the SNR daily trace
    for result in processTrace(i, elevations, S1, S2, seconds, azimuths):
      results.append(result)

  return results

def processTrace(i, elevations, S1, S2, seconds, azimuths):

  """
  def processTrace
  Processes a single trace
  """

  results = list()

  # Extract the elevation and S1 (L1) data
  for (first, last) in splitTraces(seconds):

    # Cut traces to correct indices
    traceElevation = elevations[first:last]
    traceS1 = S1[first:last]
    traceS2 = S2[first:last]
    traceAzimuth = azimuths[first:last]

    # Check if trace is OK
    if not validateTrace(traceElevation):
      printVerbose("Trace for satellite %i did not validate." % i)
      continue

    # Determine whether trace is ascending or descending
    descending = (traceElevation[0] - traceElevation[1]) > 0

    # Go over the S1 and S2 signals
    for (signal, signalType) in [(traceS1, "S1"), (traceS2, "S2")]:

      # Get the reflector height
      height = getReflectorHeight(traceElevation, signal, signalType)

      # Save the height and some auxilliary parameters for post-processing
      results.append((
        height,
        i,
        signalType,
        np.mean(traceAzimuth),
        np.mean(traceElevation),
        descending,
        len(traceElevation)
      ))

  printVerbose("SNR process completed for satellite %s with %i reflector heights." % (i, len(results)))

  return results

def getReflectorHeight(elevations, signal, type):

  """
  def getReflectorHeight
  Processes a single S1 or S2 signal
  """

  from scipy.signal import lombscargle

  # Convert dBV to Volts: unclear to what reference the dB is expressed
  # Probably does not matter as long as it is a linear scale
  # According to Shuanggen et al., 2016 this is correct:
  voltage = 10.0 ** (signal / 20.0)

  # Fit a 2nd order polynomial to the voltage data
  # Remove Pd & Pr (direct signal)
  polynomial = np.poly1d(np.polyfit(elevations, voltage, configuration.poly))
  polyvalues = polynomial(elevations)

  # Apply and subtract the fitted polynomial
  detrended = voltage - polyvalues
  
  # From Larson & Nievinski, 2012 - equation (2)
  # SNR ~ A * cos(4 * np.pi * h * (lambda ** -1) * sin(e) + phi) (eq. 2)
  # To find h, we have to identify the frequency of a function SNR = a * cos(bx + c)
  # Where:
  #    b = 4 * np.pi * h * (lambda ** -1)
  #    x = sin(e) (elevation)
  # The amplitude (A, a), phase shift (c, phi) can be ignored

  # Create sin(e) by taking the sine of the elevation
  sinElevations = np.sin(np.radians(elevations))

  # Linear space of angular frequencies to get the power at: should capture the peak
  freqs = np.linspace(1, 200, int(1E4))

  # Get the power at different periods using the Lomb Scargle algorithm for unregularly sampled data
  # Look at frequency content of SNR (dS1) as a function of sin(elevation)
  power = lombscargle(
    sinElevations,
    detrended,
    freqs,
    normalize=True
  )
  
  # Create the plot (not multiprocessing)
  if configuration.j is None and configuration.show:
    createPlot(elevations, detrended, voltage, polyvalues, freqs, power, type)

  return getHeight(freqs, power, type)


def worker(filepath):

  """
  def worker
  Worker function passed to multiprocessing unit
  """

  printVerbose("Worker %s handed file %s." % (multiprocessing.current_process(), filepath))

  # Extract the date from the filename
  date = parse(os.path.basename(filepath)[4:-10])

  # Skip all dates if requested
  if configuration.date is not None:
    if date < parse(configuration.date):
      return None

  # Extract data from the SNR file
  data = parseSNRFile(filepath, date)

  try:
    return date, processSNRFile(data)
  except Exception as e:
    print(e)
    return None


def multiprocess(files):

  """
  def multiprocess
  Enables multiple processes
  """

  # Initialize multiprocessing env
  NUMBER_OF_PROCESSES = configuration.j or multiprocessing.cpu_count()

  def initWorker():
  
    """
    def initWorker
    Initialzes a worker
    """
  
    printVerbose("Initializing new worker process: %s." % multiprocessing.current_process())
  
    # Ignore signal interrupts in the worker process
    signal.signal(signal.SIGINT, signal.SIG_IGN)

  printVerbose("Initializing pool of %i workers for processing." % NUMBER_OF_PROCESSES)

  # Create a pool (one process for each core)
  pool = multiprocessing.Pool(NUMBER_OF_PROCESSES, initWorker)

  try:
    for result in pool.imap(worker, files):
      saveResult(result)
  except KeyboardInterrupt:
    pool.terminate()

  finally:
    pool.close()
    pool.join()

def saveResult(results):

  """
  def saveResult
  Writes result of a SNR analysis to the given outfile
  """

  # No result
  if results is None:
    return

  # Unpack
  date, result = results

  printVerbose("Completed SNR process for %s with %i traces." % (date.isoformat(), len(result)))

  # Write each result to a line
  with open(configuration.output, "a") as outfile:
    for values in result:
      outfile.write("%s %.3f %i %s %.3f %.3f %i %i\n" % (date.isoformat(), *values))

def singleprocess(files):

  """
  def singleprocess
  Works through SNR files one at a time
  """

  for filepath in files:
    saveResult(worker(filepath))

def getFiles(directory):

  """
  def getFiles
  Collects files for processing from input directory
  """

  # Must be a directory
  if not os.path.isdir(directory):
    raise ValueError("The input: %s is not a directory." % directory)

  # Read the input directory and parse some metadata
  files = sorted(os.listdir(directory))
  return [os.path.join(directory, x) for x in files] 

def parseSatellites(string):

  """
  def parseSatellites
  Parses comma delimited input list to satellite range
  """

  return [int(item) for item in string.split(",")]

def parseArguments():

  """
  def parseArguments
  Parses the CMD-line arguments
  """

  import argparse

  parser = argparse.ArgumentParser(description="GNSS SNR Analysis script.")
  parser.add_argument("-j", type=int, help="Run analysis over N multiple processes. Enter 0 for the number of available cores.")
  parser.add_argument("-v", help="Run script in verbose mode.", action="store_true")
  parser.add_argument("--correction", help="Applies elevation observation angle correction.", action="store_true")
  parser.add_argument("--poly", help="Order of polynomial fit to direct power (default = 2).", default=2)
  parser.add_argument("--azimuths", type=float, help="Set azimuthal filter to accept within range (--azimuth min max).", action="append", nargs=2)
  parser.add_argument("--satellites", type=parseSatellites, help="Comma delimited list of satellites (default = all).", default=range(1, 33))
  parser.add_argument("--samples", type=int, help="Number of minimum samples for a satellite trace to be used (default = 2000).", default=2000)
  parser.add_argument("--show", help="Shows SNR plot during processing (only works without multiprocess)", action="store_true")
  parser.add_argument("--date", help="Process SNR for after single date (ISO8601).", default=None)

  # Add I/O
  parser.add_argument("input", help="Input directory of SNR files.")
  parser.add_argument("output", help="Output file to write results.")

  return parser.parse_args()

if __name__ == "__main__":

  """
  def __main__
  Extracts reflector height using GNSS reflective interferometry. 
  Example: python analyze.py -j 4 --correction --azimuths 20 70 snr outfile.dat
  """

  import sys

  # Write configuration to global scope
  configuration = parseArguments()

  printVerbose("Reading files from directory %s and writing output to %s." % (configuration.input, configuration.output))

  # Get collection of the SNR files
  files = getFiles(configuration.input) 

  printVerbose("Initializing set of %i files for processing." % len(files))

  # Check if multiprocessing is requested
  if configuration.j is not None:
    multiprocess(files)
  else:
    singleprocess(files)

  sys.exit(0)
