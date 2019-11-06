import numpy as np
import matplotlib.dates as mdates
from dateutil.parser import parse
import matplotlib.pyplot as plt

def parseData(data, type):

  dates = dict()

  for line in data:
  
    # Extact data from line
    (date, value) = line.split()
    
    # Convert to types
    if date not in dates:
      dates[date] = list()

    dates[date].append(float(value))

  pDates = []
  pVals = []

  for date in dates:

    values = np.array(dates[date])

    if len(values) > 1:
      while True:

        std = 2 * np.std(values)
        mask = np.logical_and(values < np.mean(values) + std, values > np.mean(values) - std)

        if np.all(mask):
          break

        values = values[mask]

    pDates.append(parse(date))
    pVals.append(np.mean(values))

  return pDates, pVals

if __name__ == "__main__":

  # Open the file written by analyze.py
  with open("outfile-S1.dat", "r") as infile:
    dataS1 = infile.read().split("\n")[:-1]

  # Open the file written by analyze.py
  with open("outfile-S2.dat", "r") as infile:
    dataS2 = infile.read().split("\n")[:-1]

  oneX, oneY = parseData(dataS1, "S1")
  twoX, twoY = parseData(dataS2, "S2")

  #oneY = np.convolve(oneY, np.ones((31,))/31, mode='valid')
  #twoY = np.convolve(twoY, np.ones((31,))/31, mode='valid')

  # Plot the data
  plt.scatter(oneX, oneY, c='orange')
  plt.scatter(twoX, twoY, c='blue')
  plt.show()
