"""
Example script to plot satellite
"""

import os
import datetime
import matplotlib.pyplot as plt
import numpy as np
import signal
import multiprocessing

from analyze import parseSNRFile, getPlaneCoordinates, extractMetadata

BINSIZE = (50, 50)

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

  idx = data["satellites"] == i

  azi = data["azimuths"][idx]
  eli = data["elevations"][idx]

  # No data: return an empty matrix
  if len(eli) == 0:
    return np.zeros(BINSIZE)

  x, y = getPlaneCoordinates(azi, eli)

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

if __name__ == "__main__":

  """
  def __main__
  Main function for plotting trajectories
  """

  NUMBER_OF_PROCESSES = multiprocessing.cpu_count()

  # Get collection of the SNR files
  files = sorted(os.listdir("snr"))

  # Container for all stacked heatmaps
  heatStack = np.zeros((32, *BINSIZE))

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

  for stack in heatStack:

    plt.imshow(
      stack.T,
      extent=[-1, 1, -1, 1],
      origin="lower"
    )

    plt.colorbar()
  
    plt.show()

  sys.exit(0)
