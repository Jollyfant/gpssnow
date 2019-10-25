import numpy as np
import matplotlib.dates as mdates
from dateutil.parser import parse
import matplotlib.pyplot as plt

if __name__ == "__main__":

  dates = dict()

  # Open the file written by analyze.py
  with open("outfile.dat", "r") as infile:
    data = infile.read().split("\n")[:-1]
    
  for line in data:
  
    # Extact data from line
    (date, value) = line.split()
    
    # Convert to types
    if date not in dates:
      dates[date] = list()

    dates[date].append(2.1 - float(value))

  pDates = []
  pVals = []

  for date in dates:

    if date in ["2019-03-28T00:00:00", "2019-04-17T00:00:00", "2019-04-24T09:32:00"]:
      continue

    values = np.array(dates[date])

    while True:

      std = 2 * np.std(values)
      mask = np.logical_and(values < np.mean(values) + std, values > np.mean(values) - std)

      if np.all(mask):
        break

      values = values[mask]

    pDates.append(parse(date))
    pVals.append(np.mean(values))

  plt.ylim(0)

  # Plot the data
  plt.scatter(pDates, pVals, c='orange')
  plt.show()
