#!/usr/local/bin/python3
import datetime
import os
import subprocess
import datetime

from dateutil.parser import parse

# Hoist
configuration = None

def printVerbose(message):

  """
  def printVerbose
  Prints verbose message
  """
  if configuration.v:
    print(datetime.datetime.now().isoformat(), message)

def getSP3File(date):

  """
  def getSP3File
  Convert date to spent weeks and day of week from a start date.
  """

  # Start date of GPS weeks
  GPS_START_DATE = datetime.datetime(1980, 1, 6)

  # Get the tinem difference
  delta = date - GPS_START_DATE 

  if delta.days >= 0:
    weeks = delta.days // 7
    dayOfWeek = delta.days % 7
    return os.path.join("sp3", "igs%d%d.sp3" % (weeks, dayOfWeek))
  else:
    raise ValueError("Invalid date: %s, too early." % date)

def parseArguments():

  """
  def parseArguments
  Parses the CMD-line arguments
  """

  import argparse

  # Some required arugments
  parser = argparse.ArgumentParser(description="GNSS SNR conversion script from lb2 to RINEX .o files.")
  parser.add_argument("--type", type=int, required=True, choices=[50, 66, 88, 98, 99])
  parser.add_argument("-v", help="Verbose", action="store_true")

  # Add I/O
  parser.add_argument("input", help="Input directory of raw lb2 files.")
  parser.add_argument("output", help="Output directory to write RINEX files.")

  return parser.parse_args()

if __name__ == "__main__":

  """
  def __main__
  Extracts SNR from RINEX files
  # Example python snr.py -v --type 99 RINEX-EINT tmp
  """
  # Write configuration to global scope
  configuration = parseArguments()

  if not os.path.isdir(configuration.output):
    raise ValueError("Output is not a valid directory.")
  if not os.path.isdir(configuration.input):
    raise ValueError("Input is not a valid directory.")

  for file in os.listdir(configuration.input):

    if not file.endswith(".o"):
      continue

    printVerbose("Converting RINEX file %s." % file)

    # Get the date of the file and respective SP3
    SP3File = getSP3File(parse(file[4:-6]))

    # Not available: skip
    if not os.path.exists(SP3File):
      printVerbose("Could not find SP3 file for %s." % file)
      continue

    infile = os.path.join(configuration.input, file)
    output = os.path.join(configuration.output, file + ".snr")

    subprocess.call([
      "./bin/gnssSNR",
      infile,
      output,
      SP3File,
      str(configuration.type),
    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    printVerbose("Wrote output to %s." % output)
