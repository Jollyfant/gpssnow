import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import matplotlib

from scipy.interpolate import interp1d
from dateutil.parser import parse
from dateutil.rrule import rrule, DAILY

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
  return data[abs(data - np.mean(data)) <= np.std(data)]

def runningMean(data, N):
  
  """
  def runningMean
  Simple way to calculate running mean including NaN data
  """

  # Number of samples on both sides
  TRUNC_HALF = int(0.5 * N)

  # Container for the mean values
  mean = list()

  # Go over all samples
  for i, sample in enumerate(data):

    # If the sample is NaN itself do not interpolate: keep gaps
    if np.isnan(sample):
      mean.append(np.nan)
      continue

    # Collect samples around the point but truncated edges
    minimum = max(0, i - TRUNC_HALF)
    maximum = min(len(data), i + TRUNC_HALF + 1)

    # Append the mean
    mean.append(np.nanmean(data[minimum:maximum]))

  return np.array(mean)

def parseData(data, ids):

  """
  def parseData
  Reads data from outfile
  """

  colors = list()
  dates = list()
  values = list()

  for line in data:
  
    # Extact data from line
    (date, value, id, type, azi, eli) = line.split()
    
    if int(id) != ids:
      continue

    dates.append(parse(date))
    values.append(float(value))

    # Append S2 with a different color
    if type == "S2":
      colors.append(int(id))
    else:
      colors.append(-int(id))

  dates = np.array(dates)
  values = np.array(values)
  colors = np.array(colors)

  return dates, values, colors

def plotSpatial(data):

  # Create norm azimuth between 0 and 360
  norm = plt.Normalize(0, 360)
  cmap = plt.get_cmap("hsv", 36)
  satellites = [17, 22, 28]
 
  fig, axes = plt.subplots(len(satellites))

  for i, ax in enumerate(axes):

    for year in range(2015, 2020):

      ax.axvline(x=datetime.datetime(year, 1, 1), alpha=0.25)
      ax.set_ylim(0, 4)

      # Plot satellite numbers
      text = plt.text(0.03, 0.75, satellites[i],
        horizontalalignment="center",
        color="white",
        fontsize=10,
        transform = ax.transAxes
      )

      text.set_path_effects([
        path_effects.Stroke(linewidth=4, foreground="black"),
        path_effects.Normal()
      ])

      x, y, c = extractSatellite(data, satellites[i])
      sc = ax.scatter(x, y, cmap=cmap, c=c, norm=norm)

      if i < len(satellites) - 1:
        ax.get_xaxis().set_visible(False)

  fig.text(0.075, 0.5, "Reflector Height (m)", va="center", rotation="vertical")

  # Create a colorbar for the entire plot
  cbar = fig.colorbar(sc, ax=axes, ticks=range(0, 420, 30))
  cbar.set_label("Azimuth (deg)", rotation=270, labelpad=20)

  fig.set_size_inches(10.5, 0.6 * 10.5)
  fig.savefig("color.pdf", dpi=100, bbox_inches="tight")

def extractSatellite(data, index):

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

def plotStack(data):

  """
  def plotStack
  Plot stacked satellites
  """

  # Make sure all days are included: missing days will be interpolated
  START_DATE = datetime.date(2015, 1, 1)
  END_DATE = datetime.date(2019, 5, 2)
  BASELINE = 1.9

  # Create a bucket for all dates and start appending data
  datesDict = dict()
  for dt in rrule(DAILY, dtstart=START_DATE, until=END_DATE):
    datesDict[dt.strftime("%Y-%m-%d")] = []

  for line in data:

    # Parse the data
    (date, value, id, type, azimuth, elivation) = line.split()

    # Parse as a date string
    date = parse(date).strftime("%Y-%m-%d")

    # Subtract from baseline
    value = BASELINE - float(value)

    # Add the value to the respective date bucket
    datesDict[date].append(value)

  # Reject all outliers in a single bucket
  for date in datesDict:
    datesDict[date] = rejectOutliers(np.array(datesDict[date]))

  # Extract the keys and values
  dates = np.array(list(map(lambda x: parse(x), datesDict.keys())))
  values = np.array(list(map(lambda x: np.mean(x), datesDict.values())))

  # Take a smoothed running mean for 2 weeks
  dates, values = interpolate(dates, values)
  values = runningMean(values, 9)
  values = truncate(dates, values)

  with open("snow-curve.txt", "w") as outfile:
    for date, value in zip(dates, values):
      outfile.write(date.isoformat() + " " + str(value) + "\n")
  
  # Plot and show
  plt.plot(dates, values, c="darkgray", zorder=0)
  plt.scatter(dates, values, s=4, c="gray", zorder=1)
  plt.title("Expected snow depths during winter (2015 - 2019)")
  plt.ylabel("Expected Snow Height (m)")
  plt.xlabel("Date")

  fig = plt.gcf()
  fig.set_size_inches(10.5, 0.25 * 10.5)
  fig.savefig("snow-curve.pdf", dpi=100, bbox_inches="tight")

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

    dates, values, colors = parseData(data, (i + 1))

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
  
  if PLOT_STACK:
    plotStack(data)

  if PLOT_SPATIAL:
    plotSpatial(data)
  
  if PLOT_OVERVIEW:
    plotOverview(data)
    
