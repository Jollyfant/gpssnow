"""
Example script to plot satellite x,y reflector planes
"""

import os
import datetime
import matplotlib.pyplot as plt
from matplotlib import gridspec

import sys
import numpy as np
import signal
import multiprocessing

import matplotlib.patheffects as path_effects

from analyze import parseSNRFile, getPlaneCoordinates, extractMetadata

BINSIZE = (100, 100)

def poolWorker(file):

  """
  def poolWorker
  Does work in a pool of multiprocesses
  """

  filepath, date = extractMetadata(file)

  # Open a single data file for reading
  with open(filepath, "r") as infile:
    lines = infile.read().split("\n")[:-1]

  data = parseSNRFile(lines, date)
  
  return date, np.array(list(map(lambda x: getHeatmap(data, x), range(1, 33))))

def getHeatmap(data, i):

  """ 
  def getHeatmap
  Returns stack of heatmaps for all satelittes for one file
  """ 

  # Correct satellite
  idx = data["satellites"] == i

  # Get all the azimuths and elevations
  azi = data["azimuths"][idx]
  eli = data["elevations"][idx]

  # No data: return an empty matrix
  if len(eli) == 0:
    return np.zeros(BINSIZE)

  x, y = getPlaneCoordinates(eli, azi)

  # Calculate histogram in the same [-1, 1], [-1, 1] range
  heatmap, *void = np.histogram2d(
    x, y,
    density=True,
    bins=BINSIZE,
    range=((-1, 1), (-1, 1))
  )

  return heatmap

def initWorker():

  """
  def initWorker
  Initialzes a poolWorker
  """

  print("Initializing new poolWorker process: %s." % multiprocessing.current_process())

  # Ignore signal interrupts in the poolWorker process
  signal.signal(signal.SIGINT, signal.SIG_IGN)

def plotSingle(heatStack):

  # Show single overview
  plt.imshow(
    heatStack.T, 
    extent=[-1, 1, -1, 1],
    origin="lower"
  )

  ax = plt.gca()

  text = plt.text(0.1, 0.90, 32,
    horizontalalignment="center",
    color="white",
    fontsize=10,
    transform=ax.transAxes
  )

  text.set_path_effects([
    path_effects.Stroke(linewidth=4, foreground="black"),
    path_effects.Normal()
  ])

  plt.show()

def plotOverview(heatStack):

  fig, axs = plt.subplots(8, 4)
  plt.suptitle("Reflection Tracks for GPS Satellites")

  for i, heatStack in enumerate(heatStack):

    x = i % 4
    y = int(i / 4)

    ax = axs[y, x]

    ax.imshow(
      heatStack.T,
      extent=[-1, 1, -1, 1],
      origin="lower"
    )

    ax.axis("off")

    # Plot satellite number
    text = plt.text(0.15, 0.80, str(i + 1).zfill(2),
      horizontalalignment="center",
      color="white",
      fontsize=10,
      transform = ax.transAxes
    )

    text.set_path_effects([
      path_effects.Stroke(linewidth=4, foreground="black"),
      path_effects.Normal()
    ])

  fig.set_size_inches(0.5 * 10.5, 10.5)
  fig.savefig('test2png.png', dpi=100, bbox_inches="tight")

  sys.exit(0)

if __name__ == "__main__":

  """
  def __main__
  Main function for plotting trajectories
  """

  NUMBER_OF_PROCESSES = multiprocessing.cpu_count()

  I = "snr"

  # Get collection of the SNR files
  files = sorted(os.listdir(I))

  # Container for all stacked heatmaps
  heatStack = np.zeros(BINSIZE)

  print("Initializing set of %i files for processing." % len(files))
  print("Initializing pool of %i poolWorkers for processing." % NUMBER_OF_PROCESSES)

  # Create a pool (one process for each core)
  pool = multiprocessing.Pool(NUMBER_OF_PROCESSES, initWorker)

  try:
    # Use a pool of four poolWorkers
    for date, result in pool.imap(poolWorker, files):

      print("Completed processing of file %s." % date.isoformat())

      heatStack += result

  except KeyboardInterrupt:
    pool.terminate()

  finally:
    pool.close()
    pool.join()

  #plotOverview(heatStack)
  #plotSingle(heatStack)
