import os
import sys
import subprocess

from dateutil.parser import parse

# Change this to your data directory
ROOT = "/Users/koymans/Documents/phd/data/ingv/gps"
RINEXDIR = "RINEX"

if __name__ == "__main__":

  # Make the output folder
  if not os.path.exists(RINEXDIR):
    os.makedirs(RINEXDIR)
  
  for file in sorted(os.listdir(ROOT)):
  
    if file != "ESLN201901010000.lb2":
      continue
  
    filepath = os.path.join(ROOT, file)
  
    # Write to .o and .n (nav) to RINEX folder
    with open(os.path.join(RINEXDIR, file + ".o"), "w") as outfile:
      subprocess.call([
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
        "0.0000",
        "+nav",
        os.path.join(RINEXDIR, file + ".n"),
        filepath
      ], stdout=outfile, stderr=subprocess.DEVNULL)
  