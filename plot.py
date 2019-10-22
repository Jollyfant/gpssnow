from dateutil.parser import parse
import matplotlib.pyplot as plt

if __name__ == "__main__":

  dates = []
  values = []

  # Open the file written by analyze.py
  with open("outfile.dat", "r") as infile:
    data = infile.read().split("\n")[:-1]
    
  for line in data:
  
    # Extact data from line
    (date, value) = line.split()
    
    # Convert to types
    dates.append(parse(date))
    values.append(float(value))
  
  # Plot the data
  plt.scatter(dates, values)
  plt.show()
