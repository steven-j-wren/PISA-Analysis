import numpy as np
from icecube import dataclasses, icetray
from icecube.dataio import I3File
from icecube.icetray import I3Frame, I3Units

def passed_nbi(frame):
    # TauL6 messed up for some reason. Do it manually.
    varmap = frame["TauL6_variables"]
    cc = frame["CorridorCount"].value
    z = varmap["FiniteRecoZ"]
    rho = varmap["FiniteRecoRho"]
    keepEvent = True
    if z > -225:
        keepEvent = False
    if rho > 125:
        keepEvent = False
    if z > -3 * rho:
        keepEvent = False
    keepEvent = keepEvent and (cc < 2)
    if 'Pegleg_Fit_MNHDCasc' not in frame.keys():
        keepEvent = False
    # Remove events with SRT Nchannel < 8
    hitmap = frame["SRTTWOfflinePulsesDC"].apply(frame)
    nch = 0
    for dom in hitmap.keys():
        nch += len(hitmap[dom]) > 0
    if not (nch >= 8):
        keepEvent = False

    return keepEvent
