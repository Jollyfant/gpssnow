import numpy as np
import datetime
from dateutil.parser import parse
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

# Options
PLOT_STACK = True
WRITE_OUTFILE = True

def reject_outliers(data):
  return data[abs(data - np.mean(data)) < 1 * np.std(data)]

def runningMean(N, data):
  
  """
  def runningMean
  Simple way to calculate running mean including NaN data
  """

  t = int(0.5 * N)

  mean = []

  for i in range(t, len(data) - t):
    v = []
    for j in range(i - t, i + t + 1):
      if np.isnan(data[j]):
        continue
      v.append(data[j])
    mean.append(np.sum(v) / len(v))

  return np.array(mean)

def parseData(data, ids):

  colors = list()
  dates = list()
  values = list()

  for line in data:
  
    # Extact data from line
    (date, value, id, type) = line.split()
    
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

def plotStack(data):

  """
  def plotStack
  Plot stacked satellites
  """

  dates = dict()

  for line in data:

    (date, value, id, type) = line.split()

    # Skip S2 signal
    if type == "S2":
      continue

    # Sort values by date
    if date not in dates:
      dates[date] = []

    dates[date].append(1.9 - float(value))

  for date in dates:
    dates[date] = reject_outliers(np.array(dates[date]))

  ds = list()
  vs = list()

  for date in dates:
    ds.append(parse(date))

    if(parse(date).month > 5 and parse(date).month < 11):
      vs.append(0)
    else:
      vs.append(np.mean(dates[date]))

  vs = np.array(vs)
  ds = np.array(ds)

  # Take a smoothed running mean for 1 week
  vs = runningMean(15, vs)
  ds = ds[7:-7]

  # Truncate everything below 0
  vs[vs < 0] = 0

  with open("snow-curve.txt", "w") as outfile:
    for date, value in zip(ds, vs):
      outfile.write(date.isoformat() + " " + str(value) + "\n")
  
  # Plot and show
  plt.scatter(ds, vs, s=4)
  plt.plot(ds, vs)
  plt.show()

def plotOverview(data):

  """
  def plotOverview
  Plots an overview of all satellites
  """

  fig, axs = plt.subplots(8, 4)

  for i in range(1, 33):

    dates, values, colors = parseData(data, i)

    x = (i - 1) % 4
    y = int((i - 1) / 4)

    # Plot bars for every year 
    axs[y, x].axvline(x=datetime.datetime(2015, 1, 1), alpha=0.25)
    axs[y, x].axvline(x=datetime.datetime(2016, 1, 1), alpha=0.25)
    axs[y, x].axvline(x=datetime.datetime(2017, 1, 1), alpha=0.25)
    axs[y, x].axvline(x=datetime.datetime(2018, 1, 1), alpha=0.25)
    axs[y, x].axvline(x=datetime.datetime(2019, 1, 1), alpha=0.25)

    # Plot the data
    axs[y, x].scatter(dates, values, c=colors, s=2)

    # Axes are unnecessary
    axs[y, x].get_xaxis().set_visible(False)
    axs[y, x].get_yaxis().set_visible(False)

    # Plot satellite number
    text = plt.text(0.05, 0.75, i,
      horizontalalignment='center',
      color="white",
      fontsize=12,
      transform = axs[y, x].transAxes
    )

    text.set_path_effects([
      path_effects.Stroke(linewidth=4, foreground="black"),
      path_effects.Normal()
    ])

  for ax in axs.flat:
    ax.set(xlabel="Date", ylabel="Reflector Height")

  for ax in axs.flat:
    ax.label_outer()

  plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
  plt.show()


if __name__ == "__main__":

  # Open the file written by analyze.py
  with open("outfile-good.dat", "r") as infile:
    data = infile.read().split("\n")[:-1]
  
  if PLOT_STACK:
    plotStack(data)
  else:
    plotOverview(data)
