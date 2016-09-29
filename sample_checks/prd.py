import numpy as np
from icecube import dataclasses, icetray
from icecube.dataio import I3File
from icecube.icetray import I3Frame, I3Units

def passed_prd(frame):
    # Apply final level cuts
    # First cut on PID
    keepEvent = True
    if frame["SANTA_FitType"] == 1:
        prd_pid_value = 0.6
    elif frame["SANTA_FitType"] == 2:
        prd_pid_value = 0.8
    else:
        keepEvent = False
    if frame["SANTA_PID"] >= prd_pid_value:
        keepEvent = False
    # Then cut on HLC z
    if frame["FirstHLCvertex"].pos.z > -250.0:
        keepEvent = False

    return keepEvent
