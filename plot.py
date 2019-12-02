import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import matplotlib.ticker as plticker
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
PLOT_STACK = True
PLOT_OVERVIEW = False
PLOT_SPATIAL = True
WRITE_OUTFILE = False

def rejectOutliers(data):

  """
  def rejectOutliers
  Removes values beyond 1 std in an array of floats
  """

  # No values: must be NaN
  if len(data) == 0:
    return np.nan

  # Filter everything outside 1 stdev
  return data[abs(data - np.mean(data)) <= 2 * np.std(data)]

def plotSpatial(data):

  # Create norm azimuth between 0 and 360
  norm = plt.Normalize(0, 360)
  cmap = plt.get_cmap("hsv", 36)

  # Satellites to plot
  satellites = range(0, 33)
 
  fig, axes = plt.subplots(len(satellites))

  # Go over the satellites
  for i, ax in enumerate(axes):

    ax.set_ylim(0, 4)
    ax.set_xlim(datetime.datetime(2015, 1, 1), datetime.datetime(2019, 5, 1))

    # Set the proper ticks
    ax.yaxis.set_minor_locator(plticker.MultipleLocator(base=0.5))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.xaxis.set_major_locator(mdates.YearLocator())

    # Plot satellite numbers in the corner
    text = plt.text(0.06, 0.8, "PRN " + str(satellites[i]),
      horizontalalignment="center",
      color="white",
      fontsize=12,
      transform = ax.transAxes
    )

    text.set_path_effects([
      path_effects.Stroke(linewidth=2, foreground="grey"),
      path_effects.Normal()
    ])

    # Plot the data
    x, y, c = extractSatelliteAzimuth(data, satellites[i])
    sc = ax.scatter(x, y, cmap=cmap, c=c, norm=norm, label=str(i))

    # Only bottom figure must have labels
    if i < len(satellites) - 1:
      ax.xaxis.set_major_formatter(plt.NullFormatter())
    else:
      ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

  # Text for y-axis
  fig.text(0.075, 0.5, "Reflector height (m)", va="center", rotation="vertical")

  # Create a colorbar for the entire plot
  cbar = fig.colorbar(sc, ax=axes, ticks=range(0, 420, 30))
  cbar.set_label("Average azimuth (deg)", rotation=270, labelpad=20)

  fig.set_size_inches(10.5, 0.5 * 0.33 * len(satellites) * 10.5)
  fig.savefig("color.pdf", dpi=100, bbox_inches="tight")

def extractSatelliteAzimuth(data, index):

  dates = list()
  values = list()
  azimuths = list()

  for line in data:

    # Parse the data
    (date, value, id, type, azimuth, elevation) = line.split()

    if int(id) != index:
      continue

    # Parse as a date string
    date = parse(date)
    value = float(value)
    azimuth = float(azimuth)

    dates.append(date)
    values.append(value)
    azimuths.append(azimuth)

  return (dates, values, azimuths)

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

  return meanData(N, datesDict)

def meanData(N, datesDict):

  """
  def meanData
  Takes a running mean of the cloud of scatter points
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
    
def parseData(data):

  """
  def parseData
  Plots scatter plot of all reflector heights
  """

  BASELINE = 1.9

  dates = list()
  colors = list()
  values = list()

  for line in data:

    # Parse the data and subtract baseline
    (date, value, id, type, azimuth, elivation, desc) = line.split()

    azimuth = float(azimuth)
    if azimuth < 20 or azimuth > 70:
      continue

    date = parse(date)
    value = BASELINE - float(value)

    dates.append(date)
    values.append(value)

    if type == "S1":
      colors.append(0)
      #colors.append(azimuth)
      #colors.append(int(desc))
    else:
      colors.append(1)
      #colors.append(azimuth)
      #colors.append(int(desc))

  return dates, values, colors

def plotScatter(data):

  """

  """

  (dates, values, colors) = parseData(data)

  # Take a running median!
  groupedDates, groupedValues = groupData(11, dates, values)
  groupedValues = truncate(groupedDates, groupedValues)
 
  # Convert to numpy arrays
  dates = np.array(dates)
  values = np.array(values)
  colors = np.array(colors)

  S1Dates = dates[colors == 0]
  S1Values = values[colors == 0]

  S2Dates = dates[colors == 1]
  S2Values = values[colors == 1]

  #plt.scatter(dates, values, s=1, c=colors, label="L1 frequency")
  #plt.colorbar()
  # Plot the data data
  plt.scatter(S1Dates, S1Values, s=1, c="lightgrey", label="L1 frequency")
  plt.scatter(S2Dates, S2Values, s=1, c="darkgrey", label="L2 frequency")
  plt.plot(groupedDates, groupedValues, "--", c="red", linewidth=2, label="Snow depth")

  # Layout
  plt.title("Estimated snow depth at EINT using GNSS interferometric reflectometry")
  
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

  fig = plt.gcf()
  fig.set_size_inches(10.5, 0.25 * 10.5)
  fig.savefig("snow-curve.pdf", dpi=100, bbox_inches="tight")

  with open("snow-curve.txt", "w") as outfile:
    for date, value in zip(groupedDates, groupedValues):
      outfile.write(date.isoformat() + " " + str(value) + "\n")

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

def plotOverview(data):

  """
  def plotOverview
  Plots an overview of all satellites
  """

  fig, axs = plt.subplots(2, 1)

  satellites = [19, 27]

  for _, i in enumerate(satellites):

    dates, values, colors = extractSatelliteAzimuth(data, i)

    x = _ % 4
    y = int(_ / 4)

    ax = axs[x]

    # Plot bars for every year 
    for year in range(2015, 2020):
      ax.axvline(x=datetime.datetime(year, 1, 1), alpha=0.25)

    # Plot the data
    ax.scatter(dates, values, c=colors, s=2)

    # Axes are unnecessary
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    # Plot satellite numbers
    text = plt.text(0.10, 0.75, str(i + 1).zfill(2),
      horizontalalignment="center",
      color="white",
      fontsize=10,
      transform = ax.transAxes
    )

    text.set_path_effects([
      path_effects.Stroke(linewidth=4, foreground="black"),
      path_effects.Normal()
    ])

  fig.set_size_inches(10.5, 10.5)
  fig.savefig("test2png.png", dpi=100, bbox_inches="tight")

if __name__ == "__main__":

  """
  def __main__
  Main entry point for processing SNR results to snow curve
  """

  # Open the file written by analyze.py
  with open("outfile.dat", "r") as infile:
    data = infile.read().split("\n")[:-1]
  
  plotScatter(data)
  #plotOverview(data)
  #plotSpatial(data)

  #if PLOT_SPATIAL:
    #plotSpatial(data)
  
  #if PLOT_OVERVIEW:
    #plotOverview(data)
    
