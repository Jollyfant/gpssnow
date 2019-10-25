import datetime
import os
import subprocess
from dateutil.parser import parse

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
        dayofweek = delta.days % 7
        return os.path.join("sp3", "igs%d%d.sp3" % (weeks, dayofweek))
    else:
        raise ValueError('Invalid date: %s, too early.' %date)

if __name__ == "__main__":

  GNSSSNR_OP = "99"

  # Create SNR direcrory
  if not os.path.exist("snr"):
    os.makedirs("snr")

  for file in os.listdir("RINEX"):

    # Skip .n (nav) files
    if not file.endswith(".o"):
      continue

    print("Converting RINEX file %s." % file)

    date = parse(file[4:-6])
    SP3File = getSP3File(date)

    if not os.path.exists(SP3File):
      print("Could not find SP3 file for %s." % file)
      continue

    filepath = os.path.join("RINEX", file)
    output = os.path.join("snr", file + ".snr")

    subprocess.call([
      "./bin/gnssSNR",
      filepath,
      output,
      SP3File,
      GNSSSNR_OP
    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    print("Wrote output to %s." % output)
