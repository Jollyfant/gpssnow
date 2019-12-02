#!/usr/local/bin/python3

import os
import sys
import datetime
import subprocess

from dateutil.parser import parse
from analyze import multiprocess

def printVerbose(message):

  """
  def printVerbose
  Prints verbose message
  """
  if configuration.v:
    print(datetime.datetime.now().isoformat(), message)

# Hoist configuration to global scope
configuration = None

def getSettings(type):

  """
  def getSettings
  Returns settings based on GPS type
  """

  if type == "ESLN":
    return [
      "./bin/teqc",
      "+C2",
      "-leica",
      "lb2",
      "-O.mo",
      "\"ESLN\"",
      "-O.rn",
      "355982",
      "-O.rt",
      "\"LEICA GRX1200GGPRO\"",
      "-O.an",
      "\"13286073\"",
      "-O.at",
      "\"LEIAR10\"",
      "-O.px[WGS84xyz,m]",
      "4883062.9090",
      "1306068.3008",
      "3879658.3032",
      "-O.pe[hEN,m]",
      "0.0083",
      "0.0000",
      "0.0000"
    ]

  elif type == "EINT":
    return [
      "./bin/teqc",
      "+C2",
      "-leica",
      "lb2",
      "-O.mo",
      "\"ESLN\"",
      "-O.rn",
      "455279",
      "-O.rt",
      "\"LEICA GRX1200GGPRO\"",
      "-O.an",
      "\"15072048\"",
      "-O.at",
      "\"LEIAR10\"",
      "-O.px[WGS84xyz,m]",
      "4881404.0454",
      "1307781.3919",
      "3882422.5987",
      "-O.pe[hEN,m]",
      "0.0083",
      "0.0000",
      "0.0000"
    ]

  else:
    raise ValueError("Unknown receiver type passed.")

def parseArguments():

  """
  def parseArguments
  Parses the CMD-line arguments
  """

  import argparse

  # Some required arugments
  parser = argparse.ArgumentParser(description="GNSS SNR conversion script from lb2 to RINEX .o files.")
  parser.add_argument("--type", required=True, help="Receiver type")
  parser.add_argument("-v", help="Verbose", action="store_true")

  # Add I/O
  parser.add_argument("input", help="Input directory of raw lb2 files.")
  parser.add_argument("output", help="Output directory to write RINEX files.")

  return parser.parse_args()

def worker(file):

  filepath = os.path.join(configuration.input, file)
  
  printVerbose("Converting file %s." % filepath)
 
  settings = getSettings(configuration.type) + [filepath]

  # Settings update some RINEX header fields
  with open(os.path.join(configuration.output, file + ".o"), "w") as outfile:
    subprocess.call(settings, stdout=outfile, stderr=subprocess.DEVNULL)

if __name__ == "__main__":

  """
  def __main__
  Converts raw lb2 files to RINEX .o files
  # Example python process.py -v --type EINT /Users/koymans/Documents/phd/data/ingv/int tmp
  """
  # Write configuration to global scope
  configuration = parseArguments()

  if not os.path.isdir(configuration.output):
    raise ValueError("Output is not a valid directory.") 
  if not os.path.isdir(configuration.input):
    raise ValueError("Input is not a valid directory.") 

  files = sorted(os.listdir(configuration.input))

  printVerbose("Collected %i files for processing." % len(files))

  # Check if multiprocessing is requested
  for file in files:
    worker(file)
