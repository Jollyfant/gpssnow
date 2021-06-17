import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import matplotlib.ticker as plticker
import os
import matplotlib.dates as mdates
import matplotlib
import collections

from scipy.interpolate import interp1d
from dateutil.parser import parse
from dateutil.rrule import rrule, DAILY
from itertools import chain

matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams.update({"font.size": 12})

# Options
configuration = None

def rejectOutliers(data):

  """
  def rejectOutliers
  Removes values beyond 1 std in an array of floats
  """

  # No values: must be NaN
  if len(data) == 0:
    return np.nan

  # Filter everything outside 1 stdev
  return data[abs(data - np.mean(data)) <= 1 * np.std(data)]

def plotSpatial(ids, dates, values, azimuths):

  """
  def plotSpatial
  Plots a spatial overview with azimuthal colors
  """

  # Create norm azimuth between 0 and 360
  norm = plt.Normalize(0, 360)
  cmap = plt.get_cmap("hsv", 36)

  # Satellites to plot
  fig, axes = plt.subplots(len(configuration.satellites))

  # Go over the satellites
  for i, ax in zip(configuration.satellites, axes):

    ax.set_ylim(-0.5, 1.5)
    ax.set_xlim(datetime.datetime(2015, 1, 1), datetime.datetime(2019, 5, 1))

    # Set the proper ticks
    ax.yaxis.set_minor_locator(plticker.MultipleLocator(base=0.5))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.xaxis.set_major_locator(mdates.YearLocator())

    # Plot satellite numbers in the corner
    text = plt.text(0.06, 0.8, "PRN " + str(i).zfill(2),
      horizontalalignment="center",
      color="white",
      fontsize=12,
      transform=ax.transAxes
    )

    text.set_path_effects([
      path_effects.Stroke(linewidth=2, foreground="grey"),
      path_effects.Normal()
    ])

    # Plot the data
    idx = ids == i

    x = dates[idx]
    y = values[idx]
    c = azimuths[idx]

    sc = ax.scatter(x, y, cmap=cmap, c=c, norm=norm, label=str(i))

    # Only bottom figure must have labels
    if i < len(configuration.satellites) - 1:
      ax.xaxis.set_major_formatter(plt.NullFormatter())
    else:
      ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

  # Text for y-axis
  fig.text(0.075, 0.5, "Reflector height (m)", va="center", rotation="vertical")

  # Create a colorbar for the entire plot
  cbar = fig.colorbar(sc, ax=axes, ticks=range(0, 420, 30))
  cbar.set_label("Average azimuth (deg)", rotation=270, labelpad=20)

  fig.set_size_inches(10.5, 0.5 * 0.33 * len(configuration.satellites) * 10.5)
  fig.savefig(configuration.output + "-azimuth.pdf", dpi=100, bbox_inches="tight")

def groupData(N, dates, values):

  """
  def groupData
  Create a dictionary for every day
  """

  # Begin and end date
  START_DATE = datetime.date(2015, 1, 1)
  END_DATE = datetime.date(2019, 5, 1)

  datesDict = collections.OrderedDict()

  # Sort all values per day
  for date in rrule(DAILY, dtstart=START_DATE, until=END_DATE):
    datesDict[date.strftime("%Y-%m-%d")] = list()

  for date, value in zip(dates, values):
    datesDict[date.strftime("%Y-%m-%d")].append(value)

  return medianData(N, datesDict)

def medianData(N, datesDict):

  """
  def medianData
  Takes a running mean of the cloud of curve points
  """

  meanValues = list()
  meanDates = list()

  # Number of samples on both sides
  TRUNC_HALF = int(0.5 * N)

  dates = list(datesDict.keys())
  values = list(datesDict.values())

  for i, date in enumerate(dates):

    # Minimum and maximum index of the slice (for edges)
    minimum = max(0, i - TRUNC_HALF)
    maximum = min(len(dates), i + TRUNC_HALF + 1)

    # Get the flattened slice
    valueSlice = values[minimum:maximum]
    flat = list(chain.from_iterable(valueSlice))

    # Take the median value
    if len(flat) != 0:
      value = np.median(flat)
    else:
      value = np.nan

    meanValues.append(value)
    meanDates.append(parse(date))

  return meanDates, meanValues
    
def inRange(azimuth):

  for ranges in configuration.azimuths:
    if azimuth < ranges[0] or azimuth > ranges[1]:
      return False

  return True

def parseData(data):

  """
  def parseData
  Plots curve plot of all reflector heights
  """

  dates = list()
  colors = list()
  values = list()
  azimuths = list()
  ids = list()

  for line in data:

    # Parse the data and subtract baseline
    (date, value, id, type, azimuth, elevation, desc, samples) = line.split()

    id = int(id)
    samples = int(samples)

    if samples < configuration.samples:
      continue

    # Skip some satellites
    if id not in configuration.satellites:
      continue

    azimuth = float(azimuth)

    if configuration.azimuths is not None:
      if not inRange(azimuth):
        continue

    date = parse(date)
    value = configuration.baseline - float(value)

    ids.append(id)
    dates.append(date)
    values.append(value)
    azimuths.append(azimuth)

    if type == "S1":
      colors.append(0)
    elif type == "S2":
      colors.append(1)
    else:
      raise ValueError("Unknown type.")

  return np.array(ids), np.array(dates), np.array(values), np.array(azimuths), np.array(colors)


def plotScatter(dates, values, colors):

  """
  def plotScatter
  Plots curve and fitted snow curve
  """

  # Take a running median!
  groupedDates, groupedValues = groupData(configuration.window, dates, values)

  # Truncate values
  if configuration.truncate:
    groupedValues = truncate(groupedDates, groupedValues)
 
  # Convert to numpy arrays
  dates = np.array(dates)
  values = np.array(values)
  colors = np.array(colors)

  #plt.scatter(dates, values, c=colors, s=1)
  #plt.colorbar()
  plt.show()
  #sys.exit(0)
  # Group S1 and S2 values
  S1Dates = dates[colors == 0]
  S1Values = values[colors == 0]

  S2Dates = dates[colors == 1]
  S2Values = values[colors == 1]

  ### Plot the data data
  plt.scatter(S1Dates, S1Values, s=1, c="lightgrey", label="L1 frequency")
  plt.scatter(S2Dates, S2Values, s=1, c="darkgrey", label="L2 frequency")
  plt.plot(groupedDates, groupedValues, "--", c="red", linewidth=2, label="Snow depth")

  # Layout
  plt.title("Estimated snow depth at ESLN using GNSS interferometric reflectometry")
  
  # Axes
  plt.ylabel("Snow depth (m)")
  plt.ylim(-0.5, 1.5)
  plt.xlim(groupedDates[0], groupedDates[-1])

  plt.gca().yaxis.set_minor_locator(plticker.MultipleLocator(base=0.1))
  plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
  plt.gca().xaxis.set_minor_locator(mdates.MonthLocator())
  plt.gca().xaxis.set_major_locator(mdates.YearLocator())

  # Legend
  legend = plt.gca().legend(markerscale=6, facecolor="white", loc="best", prop={"size": "x-small"}, framealpha=1)

  #plt.show()
  # Plot the snow curve
  fig = plt.gcf()
  fig.set_size_inches(10.5, 0.25 * 10.5)
  fig.savefig(configuration.output + "-curve.png", dpi=100, bbox_inches="tight")

  # Write the snow curve
  with open(configuration.output + "-curve.dat", "w") as outfile:
    for date, value in zip(groupedDates, groupedValues):
      outfile.write("%s %.3f\n" % (date.isoformat(), value))

def truncate(dates, values):

  """
  def truncate
  Truncates results to 0 sometimes
  """

  newValues = list()

  for date, value in zip(dates, values):

    # Summer months assume no snow!
    if (date.month > 5 and date.month < 11):
      value = 0

    # Cannot be below the baseline
    newValues.append(max(value, 0))

  return newValues


def interpolate(dates, values):

  """
  def interpolate
  Interpolates missing date values
  """

  # Convert dates to float timestamps
  timestamps = np.array(list(map(lambda x: x.timestamp(), dates)))

  # Filter all np.nan values: they will be interpolated
  idx = ~np.isnan(values)
  filteredTimestamps = timestamps[idx]
  values = values[idx]

  # Interpolate with a cubic spline
  spline = interp1d(
    filteredTimestamps,
    values,
    kind="cubic",
    bounds_error=False
  )

  return dates, spline(timestamps)

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

  parser = argparse.ArgumentParser(description="GNSS SNR Analysis plotting sript.")
  parser.add_argument("--truncate", help="Truncates snow depths to 0 when negative.", action="store_true")
  parser.add_argument("--curve", help="Saves curve plot of snow curve", action="store_true")
  parser.add_argument("--spatial", help="Saves spatial overview plot of satellites", action="store_true")
  parser.add_argument("--azimuths", type=float, help="Set azimuthal filter to accept within range (--azimuth min max).", action="append", nargs=2)
  parser.add_argument("--satellites", type=parseSatellites, help="Comma delimited list of satellites (default = all).", default=range(1, 33))
  parser.add_argument("--baseline", type=float, help="Baseline for reflector height", default=0)
  parser.add_argument("--samples", type=float, help="Minimum number of samples for a trace (default = 2000).", default=0)
  parser.add_argument("--window", type=float, help="Length of running median window (default = 11).", default=11)

  # Add I/O
  parser.add_argument("input", help="Input reflector height file.")
  parser.add_argument("output", help="Output name to write results.")

  return parser.parse_args()

if __name__ == "__main__":

  """
  def __main__
  Main entry point for processing SNR results to snow curve
  """

  # Write configuration to global scope
  configuration = parseArguments()

  if not os.path.isfile(configuration.input):
    raise ValueError("Input is not a valid file.")

  # Open the file written by analyze.py
  with open(configuration.input, "r") as infile:
    data = infile.read().split("\n")[:-1]
  
  (ids, dates, values, azimuths, colors) = parseData(data)

  # Make plots
  if configuration.curve:
    plotScatter(dates, values, colors)
  if configuration.spatial:
    plotSpatial(ids, dates, values, azimuths)
